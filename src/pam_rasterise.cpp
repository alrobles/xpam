// pam_rasterise.cpp -- PolygonOccSource implementation.
//
// The GDAL side rasterises feature-by-feature into a MEM dataset that
// aliases a scratch buffer sized to the geometry's pixel bounding window
// (the "MEM:::DATAPOINTER=..." pseudo-filename), so per-species work is
// O(bbox pixels) rather than O(tile pixels).  A cached per-feature
// (fid, envelope, code) table -- collected in one pass at open -- filters
// candidates per chunk; when the driver supports random reads the
// geometry is re-fetched by FID, otherwise the chunk rect is applied via
// SetSpatialFilterRect and features are re-iterated.

#include "pam_rasterise.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>
#include <ogr_spatialref.h>
#include <cpl_conv.h>
#endif

namespace pam {

#ifdef HAVE_GDAL

static void register_gdal_once() {
  static std::once_flag flag;
  std::call_once(flag, [] { GDALAllRegister(); });
}

struct PolygonOccSource::Impl {
  struct Entry {
    std::int64_t fid;
    double minx, miny, maxx, maxy;   // envelope in the LAYER's SRS
    std::int32_t code;
  };

  GDALDataset* vds = nullptr;
  OGRLayer* layer = nullptr;
  int field_idx = -1;
  bool derive = false;
  std::unordered_map<std::string, std::int32_t> code_of;
  std::vector<Entry> feats;
  std::unordered_map<std::int64_t, std::size_t> by_fid;  // filter path
  bool random_read = false;
  OGRSpatialReference* grid_srs = nullptr;   // owned
  OGRSpatialReference* layer_srs = nullptr;  // owned clone
  OGRCoordinateTransformation* to_grid = nullptr;
  OGRCoordinateTransformation* to_layer = nullptr;
  double gt[6], inv[6];
  GDALDriver* mem = nullptr;             // MEM driver handle
  std::vector<unsigned char> wbuf;

  ~Impl() {
    if (layer) layer->SetSpatialFilter(nullptr);
    if (to_grid)  OCTDestroyCoordinateTransformation(to_grid);
    if (to_layer) OCTDestroyCoordinateTransformation(to_layer);
    if (grid_srs)  grid_srs->Release();
    if (layer_srs) layer_srs->Release();
    if (vds) GDALClose(vds);
  }
};

static void apply_axis_order(OGRSpatialReference* srs) {
  if (srs) srs->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
}

PolygonOccSource::PolygonOccSource(
    const std::string& path, const std::string& layer_name,
    const std::string& field, const std::vector<std::string>& levels,
    int ncols, int nrows, const double gt[6], const std::string& srs_wkt,
    bool all_touched, int chunk_rows)
    : ncols_(ncols), nrows_(nrows),
      chunk_rows_(chunk_rows < 1 ? 256 : chunk_rows),
      all_touched_(all_touched), levels_(levels) {
  if (ncols <= 0 || nrows <= 0)
    throw std::runtime_error("PolygonOccSource: invalid grid dimensions");
  register_gdal_once();
  impl_ = new Impl();
  impl_->mem = GetGDALDriverManager()->GetDriverByName("MEM");
  if (!impl_->mem)
    throw std::runtime_error("PolygonOccSource: GDAL MEM driver missing");
  std::memcpy(impl_->gt, gt, 6 * sizeof(double));
  if (!GDALInvGeoTransform(impl_->gt, impl_->inv))
    throw std::runtime_error(
        "PolygonOccSource: env geotransform is not invertible");

  try {
    impl_->vds = (GDALDataset*)GDALOpenEx(
        path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr,
        nullptr);
    if (!impl_->vds)
      throw std::runtime_error("cannot open vector dataset: " + path);
    impl_->layer = layer_name.empty()
        ? impl_->vds->GetLayer(0)
        : impl_->vds->GetLayerByName(layer_name.c_str());
    if (!impl_->layer)
      throw std::runtime_error("layer not found in: " + path);
    impl_->field_idx =
        impl_->layer->GetLayerDefn()->GetFieldIndex(field.c_str());
    if (impl_->field_idx < 0)
      throw std::runtime_error("field '" + field + "' not found in: " + path);

    // SRS: reproject geometries only when both CRSs exist and differ.
    if (!srs_wkt.empty()) {
      impl_->grid_srs = new OGRSpatialReference();
      impl_->grid_srs->importFromWkt(srs_wkt.c_str());
      apply_axis_order(impl_->grid_srs);
    }
    const OGRSpatialReference* lref = impl_->layer->GetSpatialRef();
    if (lref) {
      impl_->layer_srs = lref->Clone();
      apply_axis_order(impl_->layer_srs);
    }
    if (impl_->grid_srs && impl_->layer_srs &&
        !impl_->layer_srs->IsSame(impl_->grid_srs)) {
      impl_->to_grid = OGRCreateCoordinateTransformation(
          impl_->layer_srs, impl_->grid_srs);
      impl_->to_layer = OGRCreateCoordinateTransformation(
          impl_->grid_srs, impl_->layer_srs);
      if (!impl_->to_grid || !impl_->to_layer)
        throw std::runtime_error(
            "cannot build coordinate transform between layer and grid SRS");
    }

    impl_->derive = levels_.empty();
    if (!impl_->derive) {
      for (std::size_t i = 0; i < levels_.size(); ++i)
        impl_->code_of.emplace(levels_[i], (std::int32_t)i);
    }

    // Pass 1: cache (fid, envelope, code) per feature -- geometry data is
    // never held here, only the envelope needed for per-chunk filtering.
    impl_->layer->ResetReading();
    bool all_fids = true;
    OGRFeature* feat;
    while ((feat = impl_->layer->GetNextFeature()) != nullptr) {
      OGRGeometry* g = feat->GetGeometryRef();
      const std::int64_t fid = feat->GetFID();
      if (fid < 0) all_fids = false;
      if (!g || !feat->IsFieldSet(impl_->field_idx)) {
        ++skipped_;
        OGRFeature::DestroyFeature(feat);
        continue;
      }
      const std::string key = feat->GetFieldAsString(impl_->field_idx);
      if (key.empty()) {
        ++skipped_;
        OGRFeature::DestroyFeature(feat);
        continue;
      }
      std::int32_t code;
      if (impl_->derive) {
        auto it = impl_->code_of.find(key);
        if (it == impl_->code_of.end()) {
          code = (std::int32_t)levels_.size();
          impl_->code_of.emplace(key, code);
          levels_.push_back(key);
        } else {
          code = it->second;
        }
      } else {
        auto it = impl_->code_of.find(key);
        if (it == impl_->code_of.end()) {
          ++skipped_;                    // field value not in `levels`
          OGRFeature::DestroyFeature(feat);
          continue;
        }
        code = it->second;
      }
      OGREnvelope e;
      g->getEnvelope(&e);
      impl_->feats.push_back({fid, e.MinX, e.MinY, e.MaxX, e.MaxY, code});
      OGRFeature::DestroyFeature(feat);
    }
    impl_->random_read = all_fids &&
        impl_->layer->TestCapability(OLCRandomRead);
    if (!impl_->random_read) {
      for (std::size_t i = 0; i < impl_->feats.size(); ++i)
        impl_->by_fid.emplace(impl_->feats[i].fid, i);
    }
  } catch (...) {
    delete impl_;
    impl_ = nullptr;
    throw;
  }
}

PolygonOccSource::~PolygonOccSource() {
  delete impl_;
}

void PolygonOccSource::fill_tile_(int y0, int h) {
  Impl* im = impl_;
  pairs_.clear();
  pos_ = 0;

  // Chunk extent in grid SRS, then in layer SRS for the feature filter.
  double xs[4], ys[4];
  const int x0s[4] = {0, ncols_, ncols_, 0};
  const int y0s[4] = {y0, y0, y0 + h, y0 + h};
  for (int k = 0; k < 4; ++k) {
    xs[k] = im->gt[0] + x0s[k] * im->gt[1] + y0s[k] * im->gt[2];
    ys[k] = im->gt[3] + x0s[k] * im->gt[4] + y0s[k] * im->gt[5];
  }
  double lx0, lx1, ly0, ly1;
  if (im->to_layer && !im->to_layer->Transform(4, xs, ys)) {
    // Rect transform failed: widen to the whole layer so nothing is missed.
    lx0 = ly0 = -std::numeric_limits<double>::max();
    lx1 = ly1 = std::numeric_limits<double>::max();
  } else {
    lx0 = lx1 = xs[0];
    ly0 = ly1 = ys[0];
    for (int k = 1; k < 4; ++k) {
      lx0 = std::min(lx0, xs[k]); lx1 = std::max(lx1, xs[k]);
      ly0 = std::min(ly0, ys[k]); ly1 = std::max(ly1, ys[k]);
    }
  }

  const auto burn_one = [&](OGRGeometry* g, std::int32_t code) {
    // Reproject to grid SRS when needed; the feature-owned geometry is
    // transformed on a clone so driver-side caches are never mutated.
    OGRGeometry* work = g;
    OGRGeometry* owned = nullptr;
    const bool linearise = !im->to_grid && g->hasCurveGeometry();
    if (im->to_grid || linearise) {
      owned = g->clone();
      if (!owned) { ++skipped_; return; }
      if (im->to_grid && owned->transform(im->to_grid) != OGRERR_NONE) {
        ++skipped_; delete owned; return;
      }
      if (owned->hasCurveGeometry()) {
        OGRGeometry* lin = owned->getLinearGeometry();
        delete owned;
        if (!lin) { ++skipped_; return; }
        owned = lin;
      }
      work = owned;
    }
    OGREnvelope e;
    work->getEnvelope(&e);                 // now in grid SRS

    // Pixel bounding window via the inverse geotransform (4 corners cover
    // the affine; rotation only widens the window conservatively).
    const double ex[4] = {e.MinX, e.MaxX, e.MaxX, e.MinX};
    const double ey[4] = {e.MinY, e.MinY, e.MaxY, e.MaxY};
    double plo = std::numeric_limits<double>::infinity();
    double phi = -plo, qlo = plo, qhi = -plo;
    for (int k = 0; k < 4; ++k) {
      const double u = im->inv[0] + im->inv[1] * ex[k] + im->inv[2] * ey[k];
      const double v = im->inv[3] + im->inv[4] * ex[k] + im->inv[5] * ey[k];
      plo = std::min(plo, u); phi = std::max(phi, u);
      qlo = std::min(qlo, v); qhi = std::max(qhi, v);
    }
    const int ix0 = std::max(0, (int)std::floor(plo));
    const int iy0 = std::max(y0, (int)std::floor(qlo));
    const int ix1 = std::min(ncols_, (int)std::ceil(phi));
    const int iy1 = std::min(y0 + h, (int)std::ceil(qhi));
    if (ix1 <= ix0 || iy1 <= iy0) { delete owned; return; }

    const int w = ix1 - ix0, wh = iy1 - iy0;
    im->wbuf.resize((std::size_t)w * wh);
    // Fresh MEM dataset over the window; MEM leaves the band uninitialised,
    // so Fill(0) before burning.  Read back via RasterIO -- the DATAPOINTER
    // alias form was removed in GDAL 3.12 for security.
    GDALDataset* m = im->mem->Create("", w, wh, 1, GDT_Byte, nullptr);
    if (!m) { delete owned; ++skipped_; return; }
    // The window's own geotransform: local pixel (0,0) is raster (ix0,iy0).
    const double wgt[6] = {im->gt[0] + ix0 * im->gt[1] + iy0 * im->gt[2],
                           im->gt[1], im->gt[2],
                           im->gt[3] + ix0 * im->gt[4] + iy0 * im->gt[5],
                           im->gt[4], im->gt[5]};
    m->SetGeoTransform(const_cast<double*>(wgt));
    GDALRasterBand* band = m->GetRasterBand(1);
    if (band->Fill(0.0) != CE_None) {
      GDALClose(m); delete owned; ++skipped_; return;
    }
    int bnum = 1;
    double burn = 1.0;
    char** opts = nullptr;
    if (all_touched_) opts = CSLAddString(opts, "ALL_TOUCHED=TRUE");
    OGRGeometryH gh = OGRGeometry::ToHandle(work);
    const CPLErr err = GDALRasterizeGeometries(
        m, 1, &bnum, 1, &gh, nullptr, nullptr, &burn, opts, nullptr, nullptr);
    CSLDestroy(opts);
    if (err != CE_None ||
        band->RasterIO(GF_Read, 0, 0, w, wh, im->wbuf.data(), w, wh,
                       GDT_Byte, 0, 0) != CE_None) {
      GDALClose(m); delete owned; ++skipped_; return;
    }
    GDALClose(m);
    delete owned;

    for (int r = 0; r < wh; ++r) {
      const std::int64_t base = (std::int64_t)(iy0 + r) * ncols_ + ix0;
      const unsigned char* row = im->wbuf.data() + (std::size_t)r * w;
      for (int c = 0; c < w; ++c)
        if (row[c])
          pairs_.push_back(((std::uint64_t)(base + c) << 32)
                           | (std::uint32_t)code);
    }
  };

  if (im->random_read) {
    for (const Impl::Entry& e : im->feats) {
      if (e.maxx < lx0 || e.minx > lx1 || e.maxy < ly0 || e.miny > ly1)
        continue;
      OGRFeature* f = im->layer->GetFeature(e.fid);
      if (!f) { ++skipped_; continue; }
      if (f->GetGeometryRef()) burn_one(f->GetGeometryRef(), e.code);
      else ++skipped_;
      OGRFeature::DestroyFeature(f);
    }
  } else {
    im->layer->SetSpatialFilterRect(lx0, ly0, lx1, ly1);
    im->layer->ResetReading();
    OGRFeature* f;
    while ((f = im->layer->GetNextFeature()) != nullptr) {
      auto it = im->by_fid.find(f->GetFID());
      if (it != im->by_fid.end() && f->GetGeometryRef())
        burn_one(f->GetGeometryRef(), im->feats[it->second].code);
      OGRFeature::DestroyFeature(f);
    }
    im->layer->SetSpatialFilter(nullptr);
  }
}

#else  // !HAVE_GDAL

void PolygonOccSource::fill_tile_(int, int) {}

#endif  // HAVE_GDAL

// ── GDAL-free bodies: the constructor is the only strict stub ─────────────

#ifndef HAVE_GDAL
PolygonOccSource::PolygonOccSource(
    const std::string&, const std::string&, const std::string&,
    const std::vector<std::string>&, int, int, const double[6],
    const std::string&, bool, int) {
  throw std::runtime_error("xpam was built without GDAL support");
}
PolygonOccSource::~PolygonOccSource() {}
#endif

std::int64_t PolygonOccSource::n_pairs() const {
  return exhausted_ ? emitted_ : -1;
}
std::int64_t PolygonOccSource::n_species() const {
  return (std::int64_t)levels_.size();
}
std::int64_t PolygonOccSource::n_skipped() const { return skipped_; }

bool PolygonOccSource::advance_() {
  while (y0_ < nrows_) {
    const int h = std::min(chunk_rows_, nrows_ - y0_);
    fill_tile_(y0_, h);
    y0_ += h;
    if (!pairs_.empty()) {
      std::sort(pairs_.begin(), pairs_.end());
      pairs_.erase(std::unique(pairs_.begin(), pairs_.end()), pairs_.end());
      pos_ = 0;
      return true;
    }
  }
  exhausted_ = true;
  return false;
}

bool PolygonOccSource::next(std::int64_t& cell, std::int32_t& sp) {
  while (pos_ >= pairs_.size()) {
    if (!advance_()) return false;
  }
  const std::uint64_t key = pairs_[pos_++];
  cell = (std::int64_t)(key >> 32);
  sp   = (std::int32_t)(key & 0xffffffffu);
  ++emitted_;
  return true;
}

} // namespace pam
