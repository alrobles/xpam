# Build a sparse environmental presence-absence matrix

Constructs the species-by-climate-bin incidence matrix (the
"environmental PAM") from a species occurrence table and a climate
table, without ever materialising a dense species-by-site or
species-by-bin matrix. Climate space is discretised onto an `nb x nb`
grid; each cell is assigned to a bin and each species is recorded as
present in the bins its cells fall in.

## Usage

``` r
pam_env(occ, env, nb = 30, xlim = NULL, ylim = NULL, max_bitset_mb = 1024)
```

## Arguments

- occ:

  A data frame with columns `cell` and `species`. `cell` must be
  matchable to `env$cell`; `species` identifies the taxon.

- env:

  A data frame with columns `cell`, `x`, `y`: the two climate
  coordinates of each cell.

- nb:

  Integer. Number of bins per climate axis. Default `30`.

- xlim, ylim:

  Numeric length-2 ranges for the binning. If `NULL`, the 0.5% and 99.5%
  quantiles of the data are used.

- max_bitset_mb:

  Numeric. Memory budget (MiB) for the bitset kernel; if the bitset
  would exceed it, the sort-based kernel is used instead. Default
  `1024`.

## Value

A list with `i`, `j` (0-based species and bin indices of the nonzero
entries), `nnz`, `n_species`, `n_bins`, `richness` (distinct species per
bin), `footprint` (bins per species), `availability` (cells per bin, the
climate-availability measure), `kernel` (the kernel actually used), and
`state_mb` (the bitset size the binning would need).

## Details

The accumulation is done in C++ with a single shared bitset
(`n_species * n_bins / 8` bytes, filled with lock-free atomic ORs), so
peak memory is independent of both the number of cells and the number of
threads, and no per-thread state has to be merged. This is what allows
the construction to run at resolutions (e.g. 1/120 degree) where the
dense path is infeasible.

The kernel is chosen from a memory budget, and the choice is reported.
If the bitset for the requested binning would exceed `max_bitset_mb`
MiB, the sort-based kernel is used instead: it has a much smaller
constant when `n_species * n_bins` is huge, but its memory is `O(P)`
(proportional to the number of (species, cell) pairs), which is what the
bitset avoids. Raise the budget, or lower `nb`, if the bitset is
affordable.

## Examples

``` r
# \donttest{
occ <- data.frame(cell = c(1, 1, 2, 3), species = c("a", "b", "a", "c"))
env <- data.frame(cell = 1:3, x = c(0.1, 0.4, 0.8), y = c(0.2, 0.5, 0.9))
pam_env(occ, env, nb = 4)
#> $n_species
#> [1] 3
#> 
#> $n_bins
#> [1] 16
#> 
#> $nnz
#> [1] 4
#> 
#> $richness
#>  [1] 2 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1
#> 
#> $footprint
#> [1] 2 1 1
#> 
#> $availability
#>  [1] 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1
#> 
#> $species
#> [1] "a" "b" "c"
#> 
#> $nb
#> [1] 4
#> 
#> $kernel
#> [1] "bitmap"
#> 
#> $state_mb
#> [1] 5.722046e-06
#> 
# }
```
