// xpam_gdal_io.h -- GDAL-backed strip reader for xpam's PamEngine.
//
// EnvReader opens D environmental rasters (either D single-band files, or one
// multi-band file whose bands are the D dimensions) and reads full-width
// strips into a pixel-major buffer: out[i * D + d] = value of pixel i on dim
// d.  Pixel-major layout feeds pam::GridBinning::assign() directly, with no
// transpose pass.
//
// The header compiles unconditionally: every method has a stub body raising
// std::runtime_error("built without GDAL support") when HAVE_GDAL is absent.
// No GDAL type appears in the public interface (the datasets live behind an
// opaque pair of pointers inside #ifdef HAVE_GDAL).
#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace pam {

class EnvReader {
public:
  // files: either one multi-band file (D = nbands) or D single-band files
  // (D = files.size()).  All bands must share raster dimensions.
  explicit EnvReader(const std::vector<std::string>& files);
  ~EnvReader();

  EnvReader(const EnvReader&)            = delete;
  EnvReader& operator=(const EnvReader&) = delete;

  int ncols() const;
  int nrows() const;
  int n_dims() const;   // D

  // Read strip rows [y0, y0 + h) x all columns into out.
  // out must hold h * ncols() * n_dims() doubles, pixel-major.
  // Band nodata values are converted to quiet NaN.
  void read_strip(int y0, int h, double* out) const;

  // Affine geotransform of the raster grid (GDAL 6-element order).
  // Throws when the dataset carries no geotransform.
  void geotransform(double gt[6]) const;

  // Raster CRS as WKT; empty string when the dataset has none.
  std::string srs_wkt() const;

  // Compile-time flag: was the package built with GDAL?
  static bool built_with_gdal() noexcept;

private:
#ifdef HAVE_GDAL
  struct Impl;
  Impl* impl_ = nullptr;
#endif
};

} // namespace pam
