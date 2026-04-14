#!/bin/bash
# Slurm job script — CDK9/CyclinT1 structure equilibration (inputs/03_equilibrate.py)
#
# Runs OpenMM minimisation + 2 ns NVT equilibration for both apo and
# holo_cyclinT1 conditions and saves 10 basis-state snapshots per condition.
#
# Estimated wall time
# -------------------
# CPU:  ~3–5 h per condition (6–10 h total)
# GPU:  ~20–40 min per condition with CUDA
#
# Usage
# -----
# From examples/openmm_cdk9_cyclinT1/:
#   sbatch scripts/hpc_equilibrate.sl
#
# Adjust the #SBATCH directives below for your cluster's partition names,
# GPU model, and account/project codes.

#SBATCH --job-name=cdk9_equilibrate
#SBATCH --output=logs/equilibrate_%j.out
#SBATCH --error=logs/equilibrate_%j.err
#SBATCH --time=12:00:00          # wall time limit (hh:mm:ss)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Uncomment and edit for GPU nodes:
# #SBATCH --gres=gpu:1
# #SBATCH --partition=gpu

# ---------------------------------------------------------------------------
# Environment setup
# ---------------------------------------------------------------------------
# Replace with your cluster's module system or conda/venv activation.
# Example for conda:
#   module load Anaconda3
#   conda activate cdk9-env
#
# Example for venv (uv):
#   source .venv/bin/activate

# Set PLATFORM to 'CUDA' if a GPU is allocated, 'CPU' otherwise.
# You can also edit PLATFORM directly in inputs/03_equilibrate.py.
export OPENMM_CPU_THREADS=${SLURM_CPUS_PER_TASK:-4}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p logs

cd "$(dirname "$0")/.."   # change to examples/openmm_cdk9_cyclinT1/

echo "Starting CDK9 equilibration on $(hostname) at $(date)"
echo "SLURM job ID: ${SLURM_JOB_ID}"

python inputs/03_equilibrate.py

echo "Equilibration complete at $(date)"
