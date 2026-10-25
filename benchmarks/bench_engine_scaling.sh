#!/bin/bash
# Driver: prep once per (side,nsp), then sweep.  Fresh process per config so
# VmHWM is a clean per-config peak RSS.
cd "$(dirname "$0")"
echo "mode,side,n_cells,n_species,threads,seconds,rss_mb,nnz"
for side in 100 316 1000 2000 4000; do
  Rscript bench_engine_scaling.R $side 500 1 prep 0.03 >&2
  for mode in engine engine_file vilela; do
    for nth in 1 8; do
      Rscript bench_engine_scaling.R $side 500 $nth $mode 0.03
    done
  done
done
for nsp in 2000 6000; do
  Rscript bench_engine_scaling.R 1000 $nsp 1 prep 0.03 >&2
  Rscript bench_engine_scaling.R 1000 $nsp 8 engine 0.03
  Rscript bench_engine_scaling.R 1000 $nsp 8 engine_file 0.03
  Rscript bench_engine_scaling.R 1000 $nsp 8 vilela 0.03
done
