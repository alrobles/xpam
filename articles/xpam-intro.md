# Getting Started with xpam

## Overview

The **xpam** package builds geographic presence-absence matrices (PAMs)
from species range polygons. A PAM is a sites-by-species binary table
where rows are grid cells and columns are species; a cell contains `1`
when the species range overlaps that cell and `0` otherwise.

PAMs are a fundamental data structure in macroecology and biogeography,
used for computing species richness, range-diversity plots,
co-occurrence analyses, and input to conservation-prioritisation
algorithms.

## Installation

``` r

# From CRAN
install.packages("xpam")

# Development version (identical code to the public repository)
# install.packages("remotes")
# remotes::install_github("alrobles/xpam")
```

## Quick Start

The core workflow has two steps:

1.  Load species range polygons as a
    [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html).
2.  Call
    [`get_geo_pam()`](https://alrobles.github.io/xpam/reference/get_geo_pam.md)
    with a desired grid resolution.

``` r

library(xpam)

# Step 1: load example data (Rodentia of Mexico)
shp <- get_rodentia_mexico()
shp

# Step 2: build the PAM at 1-degree resolution
pam_result <- get_geo_pam(shp, res = 1)
head(pam_result)
```

The result is a tibble. The first two columns (`x`, `y`) are cell
centroids; the remaining columns are species, coded as `1`/`0`.

## Adjusting Resolution

The `res` argument controls the grid cell size in the same units as the
CRS (decimal degrees for EPSG:4326): Higher resolution → more grid cells
→ finer spatial detail.

``` r

# Fine resolution (0.5 degree ~ 55 km at the equator)
pam_fine <- get_geo_pam(shp, res = 0.5)
dim(pam_fine)

# Coarse resolution (2 degrees ~ 220 km)
pam_coarse <- get_geo_pam(shp, res = 2)
dim(pam_coarse)
```

## Computing Species Richness

Once you have a PAM, species richness per cell is just a row sum:

``` r

pam_result$richness <- rowSums(pam_result[, -(1:2)])

# Quick map with base R
plot(pam_result$x, pam_result$y,
     cex  = sqrt(pam_result$richness) * 0.3,
     pch  = 15,
     col  = hcl.colors(max(pam_result$richness))[pam_result$richness],
     xlab = "Longitude",
     ylab = "Latitude",
     main = "Rodentia Species Richness — Mexico")
```

## Environmental Space

The package also provides a WorldClim PCA raster for environmental
analyses:

``` r

world_pca <- get_world_pca()
terra::plot(world_pca)
```

## The streaming engine

[`get_geo_pam()`](https://alrobles.github.io/xpam/reference/get_geo_pam.md)
builds the raster in memory, so it is meant for small inputs. For
continent- or global-scale work the package provides a streaming native
engine (C++17, OpenMP, GDAL): it rasterises geometries or bins
occurrence records straight from disk and never materialises a dense
matrix.

The engine is driven by three functions. They share the same occurrence
source (`occ`, `occ_file` or `occ_polygons`) and the same raster inputs
(`env_files`):

- [`pam_env_engine()`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)
  returns marginals only – richness per bin and footprint per species;
- [`pam_incidence()`](https://alrobles.github.io/xpam/reference/pam_incidence.md)
  returns the PAM itself as a sparse matrix;
- [`pam_jaccard()`](https://alrobles.github.io/xpam/reference/pam_jaccard.md)
  computes species-by-species (or site-by-site) Jaccard similarity from
  the same incidence.

``` r

# Environmental PAM marginals straight from raster files
eng <- pam_env_engine(
  env_files    = c("bio01.tif", "bio12.tif"),
  occ_polygons = "ranges.gpkg",
  poly_field   = "sciname",
  nbs          = c(30, 30)
)
eng$richness   # species per environmental bin
eng$footprint  # bins occupied per species
```

``` r

# The PAM itself, as a sparse units x species matrix
inc <- pam_incidence(
  env_files    = c("bio01.tif", "bio12.tif"),
  occ_polygons = "ranges.gpkg",
  poly_field   = "sciname",
  space        = "env",
  format       = "dgCMatrix"
)
dim(inc$pam)
```

``` r

# Pairwise Jaccard between species, from the same incidence
jac <- pam_jaccard(
  env_files    = c("bio01.tif", "bio12.tif"),
  occ_polygons = "ranges.gpkg",
  poly_field   = "sciname",
  space        = "env"
)
dim(jac$jaccard)
```

For geographic space pass `space = "geo"`: `env_files` then supplies
only the grid geometry and the units are `terra` cell numbers.

``` r

# Geographic PAM: raster cells as units, ranges rasterised on the fly
geo <- pam_incidence(
  env_files    = "template.tif",
  occ_polygons = "ranges.gpkg",
  poly_field   = "sciname",
  space        = "geo"
)
geo$nnz
```

The engine is built for scale, so the resource limits are explicit
rather than silent: `max_bitset_mb` caps the shared bitset, `max_out_mb`
caps the size of an incidence that may be marshalled back into R, and
`threads` controls OpenMP parallelism (`OMP_THREAD_LIMIT` is respected).

GDAL is optional.
[`has_gdal()`](https://alrobles.github.io/xpam/reference/has_gdal.md)
reports whether the installed build has it; when it does not, the engine
functions stop with an informative error and the `terra`-based dense
path keeps working.

## Using Your Own Data

You can pass any
[`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
that has a `sciname` column: each unique value of `sciname` becomes a
column in the output PAM.

``` r

library(terra)

my_shapes <- vect("path/to/my_ranges.shp")
names(my_shapes)
# [1] "sciname" "family" "order"  ...

my_pam <- get_geo_pam(my_shapes, res = 0.25)
```
