# Retrieve a WorldClim PCA Environmental Raster

Downloads a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
containing the first two principal components derived from WorldClim 2.0
bioclimatic variables at global extent.

## Usage

``` r
get_world_pca()
```

## Source

WorldClim 2.0 project: <https://www.worldclim.org/>

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object, or `"Can't access data."` if the remote server is unreachable.

## Examples

``` r
# \donttest{
ras <- get_world_pca()
# Guarded: the function returns a character string when the data
# server is unreachable, so the example also runs offline
if (!is.character(ras)) {
  print(ras)
  terra::plot(ras)
}
#> class       : SpatRaster
#> size        : 1042, 2160, 2  (nrow, ncol, nlyr)
#> resolution  : 0.1666667, 0.1666667  (x, y)
#> extent      : -180, 180, -90, 83.66667  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84
#> source(s)   : memory
#> names       :        PC1,       PC2
#> min values  : -23.130791, -4.089408
#> max values  :   5.695019,   18.6616

# }
```
