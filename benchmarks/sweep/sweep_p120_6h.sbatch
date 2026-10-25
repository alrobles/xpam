#!/usr/bin/env bash
#SBATCH --job-name=xpam_sweep_p120_6h
#SBATCH --partition=sixhour
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=05:59:00
#SBATCH --output=/beegfs/a474r867/xpam_sweep/logs/sweep_p120_6h_%j.out
#SBATCH --error=/beegfs/a474r867/xpam_sweep/logs/sweep_p120_6h_%j.err

# Opportunistic p120 run inside the 6-hour partition while the full
# sweep waits on kbs (job 31837115).  Only the informative configs:
# poly threads=1 (writes the pair table), file x3, mem x1 (reports
# infeasible by design via the 24 GB gate).  Polygon rasterisation
# measured thread-flat at t10, so poly 8/24 add nothing.  A separate
# XPOC name avoids any collision with the kbs job's occ_p120.xpoc.
export RES_LIST="p120"
export ONLY="p120:poly:1,p120:file:*,p120:mem:1"
export XPOC=/beegfs/a474r867/xpam_sweep/work/occ_p120_s6.xpoc
bash /beegfs/a474r867/xpam_sweep/sweep.sh
