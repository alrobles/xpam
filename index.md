# xpam

**xpam** builds presence-absence matrices (PAMs) from species range
polygons or occurrence records. A PAM is a sites-by-species binary
table: rows are sites (raster cells or environmental bins), columns are
species, and each entry is `1` (present) or `0` (absent). PAMs are a
foundational data structure in macroecology and biogeography for
computing species richness, range-diversity plots, co-occurrence
analysis, and conservation prioritisation.

The package offers two scales of operation:

- a **dense small-input path**
  ([`get_geo_pam()`](https://alrobles.github.io/xpam/reference/get_geo_pam.md))
  that rasterises a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  onto a grid and returns a tibble;
- a **streaming native engine** (C++17, OpenMP, GDAL optional) that
  never materialises a dense matrix —
  [`pam_incidence()`](https://alrobles.github.io/xpam/reference/pam_incidence.md)
  returns the sparse PAM,
  [`pam_env_engine()`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)
  returns marginals, and
  [`pam_jaccard()`](https://alrobles.github.io/xpam/reference/pam_jaccard.md)
  computes species- or site-wise Jaccard similarity directly from the
  same incidence. Proven at 22.5e9 occurrence pairs.

## Installation

``` r

# From CRAN
install.packages("xpam")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("alrobles/xpam")
```

GDAL is optional. When present at build time
([`has_gdal()`](https://alrobles.github.io/xpam/reference/has_gdal.md)
reports it) the native raster engine links and streams vector sources
directly; without GDAL the package still installs and the terra-based
path works.

## Quick Start

``` r

library(xpam)

# Load example data: Rodentia of Mexico (downloaded on demand)
shp <- get_rodentia_mexico()

# Build a PAM at 1-degree resolution
pam_result <- get_geo_pam(shp, res = 1)
head(pam_result)
#> # A tibble: 6 × 42
#>       x     y Chaetodipus arenarius ...
#>   <dbl> <dbl>                 <dbl>
#> 1 -118.  23.5                     0
#> ... 
```

A fully local example that needs no download:

``` r

square <- function(x0, y0, size = 2) {
  m <- matrix(c(x0, y0, x0 + size, y0,
                x0 + size, y0 + size, x0, y0 + size),
              ncol = 2, byrow = TRUE)
  terra::vect(m, type = "polygons", crs = "EPSG:4326")
}
ranges <- rbind(square(0, 0), square(2, 2))
ranges$sciname <- c("species_a", "species_b")
get_geo_pam(ranges, res = 1)
```

## Main Functions

| Function | Description |
|----|----|
| [`get_geo_pam()`](https://alrobles.github.io/xpam/reference/get_geo_pam.md) | Dense PAM from a `SpatVector` at a given resolution |
| [`pam_incidence()`](https://alrobles.github.io/xpam/reference/pam_incidence.md) | Sparse PAM (`Matrix::dgCMatrix` or `i`/`j` triplets) from any occurrence source |
| [`pam_env_engine()`](https://alrobles.github.io/xpam/reference/pam_env_engine.md) | Streaming environmental PAM marginals from raster files |
| [`pam_env()`](https://alrobles.github.io/xpam/reference/pam_env.md) | Dense small-input environmental PAM |
| [`pam_jaccard()`](https://alrobles.github.io/xpam/reference/pam_jaccard.md) | Species- or site-wise Jaccard similarity (optional thresholded edge list) |
| [`has_gdal()`](https://alrobles.github.io/xpam/reference/has_gdal.md) | Report whether the build has GDAL support |
| [`get_rodentia_mexico()`](https://alrobles.github.io/xpam/reference/get_rodentia_mexico.md) | Example data: rodent ranges for Mexico |
| [`get_world_pca()`](https://alrobles.github.io/xpam/reference/get_world_pca.md) | WorldClim bioclimatic PCA raster (2 layers) |

## Documentation

Full documentation, reference and the introductory article live at
<https://alrobles.github.io/xpam/>. Vignettes are published there and
are deliberately not shipped in the CRAN tarball.

## Citation

``` r

citation("xpam")
```

## License

GPL (\>= 3)
