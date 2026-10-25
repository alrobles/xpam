# Check whether the package was built with GDAL support

Check whether the package was built with GDAL support

## Usage

``` r
has_gdal()
```

## Value

Logical scalar: `TRUE` when the native GDAL engine
([`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md))
is available.

## Examples

``` r
has_gdal()
#> [1] TRUE
```
