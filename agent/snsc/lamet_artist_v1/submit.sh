#!/bin/bash
#SBATCH -J lamet_meff
#SBATCH -p cpu6248R,cpueicc,i72c512g
#SBATCH --output=/public/home/zhangxin/lattice-pdf/agent/snsc/logs/lamet_%j.out
#SBATCH --error=/public/home/zhangxin/lattice-pdf/agent/snsc/logs/lamet_%j.err
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --chdir=/public/home/zhangxin/lattice-pdf/agent/snsc

set -euo pipefail
source /public/home/zhangxin/miniconda3/etc/profile.d/conda.sh && conda activate zhangxin-snsc
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

echo "=== lamet-agent Correlator Artist ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

python -u lamet_artist_v1/meff_analysis.py \
    --data-dir artifacts \
    --conf-id 46000 \
    --Pz-list="-2,-3,-4,-5,-6" \
    --nstate 2 \
    --output-dir meff_analysis 2>&1

echo "Finished: $(date)"
