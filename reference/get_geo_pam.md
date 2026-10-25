# Build a Presence-Absence Matrix from Species Range Polygons

Rasterizes a set of species range polygons onto a regular grid and
returns a sites-by-species presence-absence matrix (PAM) as a tibble.

## Usage

``` r
get_geo_pam(shp, res = 0.5, ...)
```

## Arguments

- shp:

  A
  [`SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  with a `sciname` attribute column identifying each species.

- res:

  Numeric. Grid cell size in the same units as the coordinate reference
  system of `shp` (typically decimal degrees). Default `0.5`.

- ...:

  Additional arguments passed to
  [`rasterize`](https://rspatial.github.io/terra/reference/rasterize.html).

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns `x`, `y` (cell centroids) followed by one column per species
containing `1` (present) or `0` (absent).

## Examples

``` r
# \donttest{
# A small, fully local example: two square ranges, no download needed
square <- function(x0, y0, size = 2) {
  m <- matrix(c(x0, y0, x0 + size, y0,
                x0 + size, y0 + size, x0, y0 + size),
              ncol = 2, byrow = TRUE)
  terra::vect(m, type = "polygons", crs = "EPSG:4326")
}
ranges <- rbind(square(0, 0), square(2, 2))
ranges$sciname <- c("species_a", "species_b")
pam <- get_geo_pam(ranges, res = 1)
head(pam)
#> # A tibble: 6 × 4
#>       x     y species_a species_b
#>   <dbl> <dbl>     <dbl>     <dbl>
#> 1   0.5   0.5         1         0
#> 2   0.5   1.5         1         0
#> 3   1.5   0.5         1         0
#> 4   1.5   1.5         1         0
#> 5   2.5   2.5         0         1
#> 6   2.5   3.5         0         1

# The packaged example data are downloaded on demand, so guard for
# machines without network access (the function returns a character
# string instead of a SpatVector when the server is unreachable)
shp <- get_rodentia_mexico()
if (!is.character(shp)) {
  pam_mex <- get_geo_pam(shp, res = 1)
  head(pam_mex)
}
#> # A tibble: 6 × 27
#>       x     y `Orthogeomys grandis` `Peromyscus gambelii` Reithrodontomys fulv…¹
#>   <dbl> <dbl>                 <dbl>                 <dbl>                  <dbl>
#> 1 -117.  32.0                     0                     0                      0
#> 2 -116.  30.0                     0                     0                      0
#> 3 -116.  31.0                     0                     0                      0
#> 4 -116.  32.0                     0                     0                      0
#> 5 -115.  30.0                     0                     0                      0
#> 6 -115.  32.0                     0                     0                      0
#> # ℹ abbreviated name: ¹​`Reithrodontomys fulvescens`
#> # ℹ 22 more variables: `Reithrodontomys megalotis` <dbl>,
#> #   `Sciurus alleni` <dbl>, `Sciurus niger` <dbl>, `Sigmodon alleni` <dbl>,
#> #   `Sigmodon arizonae` <dbl>, `Sigmodon hirsutus` <dbl>,
#> #   `Sigmodon hispidus` <dbl>, `Sigmodon leucotis` <dbl>,
#> #   `Sigmodon mascotensis` <dbl>, `Sigmodon ochrognathus` <dbl>,
#> #   `Sigmodon planifrons` <dbl>, `Sigmodon toltecus` <dbl>, …
# }
```
