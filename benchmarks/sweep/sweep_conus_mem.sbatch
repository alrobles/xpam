#!/usr/bin/env bash
#SBATCH --job-name=xpam_conus_mem
#SBATCH --partition=sixhour
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/beegfs/a474r867/xpam_sweep/logs/sweep_conus_mem_%j.out
#SBATCH --error=/beegfs/a474r867/xpam_sweep/logs/sweep_conus_mem_%j.err

# Records the conus `mem` row: the 24 GB gate in sweep_run.R reads only
# the XPOC header, reports the resident-pair-table config infeasible
# (~37 GB at 9.2e8 pairs) and prints the NA row -- the designed result.
export RES_LIST="conus"
export ONLY="conus:mem:*"
bash /beegfs/a474r867/xpam_sweep/sweep.sh
