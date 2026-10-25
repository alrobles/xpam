#!/usr/bin/env bash
#SBATCH --job-name=xpam_sweep_conus
#SBATCH --partition=sixhour
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=05:59:00
#SBATCH --output=/beegfs/a474r867/xpam_sweep/logs/sweep_conus_%j.out
#SBATCH --error=/beegfs/a474r867/xpam_sweep/logs/sweep_conus_%j.err

# Conterminous-USA window at 1/120 deg: same resolution and code path
# as the global p120 run but ~2.1e7 cells, so the full config matrix
# fits one six-hour slot while the global job waits on kbs.
export RES_LIST="conus"
bash /beegfs/a474r867/xpam_sweep/sweep.sh
