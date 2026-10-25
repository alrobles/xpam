#!/usr/bin/env bash
#SBATCH --job-name=xpam_sweep_t10
#SBATCH --partition=sixhour
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=05:00:00
#SBATCH --output=/beegfs/a474r867/xpam_sweep/logs/sweep_t10_%j.out
#SBATCH --error=/beegfs/a474r867/xpam_sweep/logs/sweep_t10_%j.err

export RES_LIST="t10"
bash /beegfs/a474r867/xpam_sweep/sweep.sh
