#!/usr/bin/env bash
#SBATCH --job-name=xpam_sweep_p120
#SBATCH --partition=kbs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=96:00:00
#SBATCH --output=/beegfs/a474r867/xpam_sweep/logs/sweep_p120_%j.out
#SBATCH --error=/beegfs/a474r867/xpam_sweep/logs/sweep_p120_%j.err

# Global 1/120-deg sweep on the long partition.  Polygon rasterisation
# measured thread-flat at t10 (1371/1388/1358 s) and conus
# (4646/4572/4544 s), so poly runs ONCE at threads=1 -- that run also
# writes occ_p120.xpoc, which the three `file` runs then stream.
# mem:1 reports infeasible by design via the 24 GB gate (the resident
# pair table would be ~100+ GB).  Time 48h -> 96h: conus (2.1e7 cells,
# ~9.2e8 pairs) took ~76 min of rasterisation; global scales with
# polygon coverage and may take 1-2 days.  Previous attempts: sixhour
# job 31849178 hit the 6h limit; kbs job 31837115 ran sweep.sh mid-edit
# and produced an empty CSV.
export RES_LIST="p120"
export ONLY="p120:poly:1,p120:file:*,p120:mem:1"
bash /beegfs/a474r867/xpam_sweep/sweep.sh
