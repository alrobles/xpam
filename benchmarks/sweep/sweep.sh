#!/usr/bin/env bash
# sweep.sh -- prep + sequential config grid for the xpam real-data sweep.
# Run on the cluster under SLURM (see sweep.sbatch); every config is a
# fresh Rscript process so /proc/self/status VmHWM is a clean peak RSS.

set -euo pipefail
ROOT=/beegfs/a474r867/xpam_sweep
WORK=$ROOT/work
RES_OUT=$ROOT/results
YEAR=${YEAR:-2020}
IUCN=/home/a474r867/work/gbifdata/iucn_shapefile

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xpam
mkdir -p "$WORK" "$RES_OUT" "$ROOT/logs"

CSV=$RES_OUT/sweep_${YEAR}_$(date +%Y%m%d_%H%M%S).csv
echo "mode,res,n_cells,n_species,threads,seconds,rss_mb,nnz,n_skipped" > "$CSV"

# ── prep: merge the two IUCN mammal parts into one GeoPackage ───────────
# IUCN protocol: presence = 1 (extant) and origin = 1 (native) only.
if [[ ! -s $WORK/mammals.gpkg ]]; then
  ogr2ogr -f GPKG "$WORK/mammals.gpkg" "$IUCN/MAMMALS_PART1.shp" \
          -nln mammals -nlt PROMOTE_TO_MULTI \
          -where "presence = 1 AND origin = 1"
  ogr2ogr -f GPKG -update -append "$WORK/mammals.gpkg" \
          "$IUCN/MAMMALS_PART2.shp" -nln mammals \
          -where "presence = 1 AND origin = 1"
  ogrinfo "$WORK/mammals.gpkg" -dialect SQLITE \
          -sql "SELECT COUNT(*) n_features, COUNT(DISTINCT sci_name) n_species FROM mammals" \
          >&2
fi

# RES_LIST selects which resolutions this invocation runs (default both).
RES_LIST=${RES_LIST:-"t10 p120"}

# ── prep: per-resolution bin bounds (0.5--99.5% band quantiles) + ncells ─
# NOTE: terra::global(quantile) materialises every cell value -> OOM at
# 1/120 deg under a SLURM cgroup.  A regular spatial sample of <=2e6
# cells gives the same bounds at bounded memory.
for R in $RES_LIST; do
  if [[ ! -s $WORK/meta_$R.rds ]]; then
    Rscript - "$R" "$YEAR" <<'EOF'
suppressMessages(library(terra))
args <- commandArgs(trailingOnly = TRUE); rid <- args[1]; yr <- args[2]
DIV   <- "/beegfs/a474r867/diversity_work"
STAGE <- "/beegfs/a474r867/stage1_global_hybrid/bioclim_separate_1km_global_1986-2025"
ef <- switch(rid,
  t10  = file.path(DIV, "coarse50", sprintf("bio%02d_%s.tif", c(1,12), yr)),
  p120 = file.path(STAGE, yr, sprintf("bio%02d_%s.tif", c(1,12), yr)),
  conus = file.path("/beegfs/a474r867/xpam_sweep/work",
                    sprintf("conus_bio%02d_%s.tif", c(1,12), yr)))
r <- rast(ef)
s <- spatSample(r, size = min(2e6, ncell(r)), method = "regular",
                na.rm = TRUE, warn = FALSE)
qs <- apply(as.matrix(s), 2L, stats::quantile,
            probs = c(0.005, 0.995), na.rm = TRUE)
saveRDS(list(mins = qs[1, ], maxs = qs[2, ],
             ncells = ncell(r), ncols = ncol(r), nrows = nrow(r)),
      file.path("/beegfs/a474r867/xpam_sweep/work",
                sprintf("meta_%s.rds", rid)))
EOF
  fi
done

# ── config grid: poly first (writes the XPOC the other modes stream) ────
# ONLY (comma-separated res:mode:threads globs) filters configs;
# a failing config logs and does not abort the sweep.
ONLY=${ONLY:-"*"}
run() { # <res> <mode> <threads>
  local key="$1:$2:$3" match=0 p
  local pats=()
  IFS=',' read -ra pats <<< "$ONLY"   # split without pathname expansion
  for p in "${pats[@]}"; do
    [[ "$key" == $p ]] && { match=1; break; }
  done
  (( match )) || return 0
  echo ">>> $1 $2 threads=$3" >&2
  Rscript "$ROOT/sweep_run.R" "$1" "$2" "$3" "$YEAR" >> "$CSV" ||
    echo "!!! FAILED: $key" >&2
}

if [[ " $RES_LIST " == *" t10 "* ]]; then
  # 1/10 deg -- all three modes incl. in-memory
  for th in 1 8 24; do run t10 poly $th; done
  for th in 1 8 24; do run t10 file $th; done
  for th in 1 8 24; do run t10 mem  $th; done
fi

if [[ " $RES_LIST " == *" p120 "* ]]; then
  # 1/120 deg -- polygon rasterise once per thread count, then file/mem
  for th in 1 8 24; do run p120 poly $th; done
  for th in 1 8 24; do run p120 file $th; done
  run p120 mem 1   # pair table at 1/120 may be large; one in-memory check
fi

if [[ " $RES_LIST " == *" conus "* ]]; then
  # conterminous USA window at 1/120 deg (~2.1e7 cells) -- same code,
  # regional extent so poly fits a six-hour slot
  for th in 1 8 24; do run conus poly $th; done
  for th in 1 8 24; do run conus file $th; done
  run conus mem 1
fi

echo "done -> $CSV" >&2
cat "$CSV" >&2
