// pam_rasterise.h -- streaming polygon -> (species, cell) rasteriser
// (Phase D).
//
// PolygonOccSource is an OccSource whose pairs come from
// GDALRasterizeGeometries applied to row-chunks of the target grid.  A
// first pass caches (fid, envelope, species code) per feature, so memory
// is O(n_features), not O(geometry).  Per chunk, candidate geometries are
// burned one by one into a scratch MEM band sized to each polygon's own
// bounding window (never a full-tile buffer), and the emitted cells are
// sorted inside the chunk.  Chunks advance in row order, so the emitted
// stream is globally cell-sorted -- the contract the engine's merge pass
// relies on.
//
// Species codes are the 0-based index into `levels`; an empty `levels`
// derives species order from first occurrence in the layer.  When the
// vector SRS differs from the grid SRS, geometries are reprojected via an
// OGRCoordinateTransformation; when either SRS is missing they are assumed
// identical.  Dumping the emitted stream to an XPOC file lives in
// DumpOccSource (applied post-merge), not here.
//
// Compiles unconditionally; without HAVE_GDAL the constructor throws.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "pam_occ_io.h"

namespace pam {

class PolygonOccSource : public OccSource {
public:
  // path:   any OGR-readable vector dataset (shp, gpkg, ...).
  // layer:  layer name; "" selects the first layer.
  // field:  attribute identifying the species; values are read as strings
  //         and matched against `levels`.
  // levels: species order -> code (index); empty derives it in file order.
  // ncols/nrows/gt: target raster grid (usually the env raster's).
  // srs_wkt: grid CRS; "" assumes the layer's own coordinates.
  // all_touched: GDAL ALL_TOUCHED burn rule (vs pixel-centre).
  // chunk_rows: rows per rasterise/emit chunk.
  PolygonOccSource(const std::string& path, const std::string& layer,
                   const std::string& field,
                   const std::vector<std::string>& levels,
                   int ncols, int nrows, const double gt[6],
                   const std::string& srs_wkt, bool all_touched,
                   int chunk_rows);
  ~PolygonOccSource() override;

  PolygonOccSource(const PolygonOccSource&)            = delete;
  PolygonOccSource& operator=(const PolygonOccSource&) = delete;

  // Total pairs is not known until the stream is exhausted (returned -1
  // while still draining).
  std::int64_t n_pairs()   const override;
  std::int64_t n_species() const override;
  std::int64_t n_skipped() const override;

  // Species levels actually used (derived or as given).
  const std::vector<std::string>& species_levels() const { return levels_; }

  bool next(std::int64_t& cell, std::int32_t& sp) override;

private:
  bool advance_();          // fill pairs_ from the next row-chunk
  void fill_tile_(int y0, int h);

  int ncols_, nrows_, chunk_rows_;
  int y0_ = 0;
  bool all_touched_;
  bool exhausted_ = false;
  std::vector<std::string> levels_;
  std::vector<std::uint64_t> pairs_;   // packed (cell << 32) | sp, sorted
  std::size_t pos_ = 0;
  std::int64_t emitted_ = 0, skipped_ = 0;
#ifdef HAVE_GDAL
  struct Impl;
  Impl* impl_ = nullptr;
#endif
};

} // namespace pam
