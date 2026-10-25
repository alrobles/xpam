// xpam_gdal_io.cpp -- EnvReader implementation.
// Bodies guarded by HAVE_GDAL; the #else side of every method throws a
// call-time error so the package always compiles and loads.

#include "xpam_gdal_io.h"

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <cmath>
#include <limits>
#include <mutex>
#endif

namespace pam {

#ifdef HAVE_GDAL

struct EnvReader::Impl {
  // (dataset, band_index) per dimension: one entry per dim d, so a multi-band
  // file and a list of single-band files share one code path.
  std::vector<std::pair<GDALDataset*, int>> bands;
  std::vector<GDALDataset*> owned;                 // distinct handles to close
  std::vector<double> nodata;                      // per dim; NaN if unset
  std::vector<bool>   has_nodata;
  double gt[6] = {0, 1, 0, 0, 0, 1};               // identity fallback
  bool has_gt = false;
  std::string srs;
  int nc = 0, nr = 0;

  ~Impl() { for (GDALDataset* d : owned) if (d) GDALClose(d); }
};

static void register_gdal_once() {
  static std::once_flag flag;
  std::call_once(flag, [] { GDALAllRegister(); });
}

EnvReader::EnvReader(const std::vector<std::string>& files) {
  register_gdal_once();
  impl_ = new Impl();
  try {
    std::vector<std::pair<std::string, int>> spec;   // (path, band)
    if (files.size() == 1) {
      GDALDataset* ds = (GDALDataset*)GDALOpen(files[0].c_str(), GA_ReadOnly);
      if (!ds) throw std::runtime_error("EnvReader: cannot open '" + files[0] + "'");
      impl_->owned.push_back(ds);
      const int nb = ds->GetRasterCount();
      if (nb < 1) throw std::runtime_error("EnvReader: '" + files[0] + "' has no bands");
      for (int b = 1; b <= nb; ++b) spec.push_back({files[0], b});
    } else {
      for (const auto& f : files) spec.push_back({f, 1});
    }
    for (const auto& s : spec) {
      GDALDataset* ds;
      if (s.first == files[0] && files.size() == 1) {
        ds = impl_->owned[0];
      } else {
        ds = (GDALDataset*)GDALOpen(s.first.c_str(), GA_ReadOnly);
        if (!ds) throw std::runtime_error("EnvReader: cannot open '" + s.first + "'");
        impl_->owned.push_back(ds);
      }
      if (impl_->bands.empty()) {
        impl_->nc = ds->GetRasterXSize();
        impl_->nr = ds->GetRasterYSize();
        impl_->has_gt = ds->GetGeoTransform(impl_->gt) == CE_None;
        const char* proj = ds->GetProjectionRef();
        if (proj) impl_->srs = proj;
      } else if (ds->GetRasterXSize() != impl_->nc ||
                 ds->GetRasterYSize() != impl_->nr) {
        throw std::runtime_error(
          "EnvReader: raster dimensions differ across env files");
      }
      impl_->bands.push_back({ds, s.second});
      GDALRasterBand* band = ds->GetRasterBand(s.second);
      int has_nd = 0;
      double nd = band->GetNoDataValue(&has_nd);
      impl_->has_nodata.push_back(has_nd != 0);
      impl_->nodata.push_back(nd);
    }
  } catch (...) {
    delete impl_; impl_ = nullptr;
    throw;
  }
}

EnvReader::~EnvReader() { delete impl_; }

int EnvReader::ncols() const { return impl_->nc; }
int EnvReader::nrows() const { return impl_->nr; }
int EnvReader::n_dims() const { return (int)impl_->bands.size(); }

void EnvReader::read_strip(int y0, int h, double* out) const {
  const int nc = impl_->nc, D = (int)impl_->bands.size();
  const std::size_t npix = (std::size_t)h * nc;
  for (int d = 0; d < D; ++d) {
    GDALRasterBand* band = impl_->bands[d].first->GetRasterBand(impl_->bands[d].second);
    // Write band d straight into its pixel-major slot.  RasterIO strides are
    // in BYTES: pixel stride D * sizeof(double), line stride D * nc *
    // sizeof(double) -- no transpose pass needed.
    CPLErr err = band->RasterIO(GF_Read, 0, y0, nc, h,
                                out + d, nc, h, GDT_Float64,
                                (GSpacing)(D * sizeof(double)),
                                (GSpacing)(D * nc * sizeof(double)),
                                nullptr);
    if (err != CE_None)
      throw std::runtime_error("EnvReader::read_strip: RasterIO failed");
    if (impl_->has_nodata[d]) {
      const double nd = impl_->nodata[d];
      for (std::size_t i = 0; i < npix; ++i) {
        double& v = out[i * D + d];
        if (v == nd) v = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void EnvReader::geotransform(double gt[6]) const {
  if (!impl_->has_gt)
    throw std::runtime_error("EnvReader: env raster has no geotransform");
  for (int i = 0; i < 6; ++i) gt[i] = impl_->gt[i];
}

std::string EnvReader::srs_wkt() const { return impl_->srs; }

bool EnvReader::built_with_gdal() noexcept { return true; }

#else  // !HAVE_GDAL -- stubs: compile always, fail at call time

EnvReader::EnvReader(const std::vector<std::string>&) {
  throw std::runtime_error("xpam was built without GDAL support");
}
EnvReader::~EnvReader() {}
int EnvReader::ncols() const { throw std::runtime_error("xpam was built without GDAL support"); }
int EnvReader::nrows() const { throw std::runtime_error("xpam was built without GDAL support"); }
int EnvReader::n_dims() const { throw std::runtime_error("xpam was built without GDAL support"); }
void EnvReader::read_strip(int, int, double*) const {
  throw std::runtime_error("xpam was built without GDAL support");
}
void EnvReader::geotransform(double*) const {
  throw std::runtime_error("xpam was built without GDAL support");
}
std::string EnvReader::srs_wkt() const {
  throw std::runtime_error("xpam was built without GDAL support");
}
bool EnvReader::built_with_gdal() noexcept { return false; }

#endif

} // namespace pam
