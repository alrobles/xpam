# Retrieve Rodentia Range Maps for Mexico

Downloads a pre-built
[`SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
containing IUCN range polygons for all Rodentia species recorded in
Mexico.

## Usage

``` r
get_rodentia_mexico()
```

## Source

Robles-Fernandez, A. L., Lira-Noriega, A., & Martinez-Meyer, E. (2023).
A phylogeny-informed characterisation of global tetrapod traits
addresses data gaps and biases.
[doi:10.1101/2023.03.04.531098](https://doi.org/10.1101/2023.03.04.531098)

## Value

A
[`SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
object, or `"Can't access data."` if the remote server is unreachable.

## Examples

``` r
# \donttest{
shp <- get_rodentia_mexico()
# Guarded: the function returns a character string when the data
# server is unreachable, so the example also runs offline
if (!is.character(shp)) {
  print(shp)
}
#> class       : SpatVector
#> geometry    : polygons
#> dimensions  : 247, 1  (geometries, attributes)
#> extent      : -117.1277, -86.81307, 14.54421, 32.72083  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> names       :                    sciname
#> type        :                      <chr>
#> values      :  Ammospermophilus harrisii
#>               Ammospermophilus interpres
#>                Ammospermophilus leucurus
#>               ...
# }
```
