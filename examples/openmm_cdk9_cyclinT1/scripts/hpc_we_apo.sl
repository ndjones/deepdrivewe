#!/bin/bash
# Slurm job script — CDK9 apo weighted ensemble (Phase 1)
#
# Runs the CDK9 apo weighted ensemble using Academy agents.
# Each CDK9SimulationAgent worker occupies one GPU.  Set num_workers in
# config_apo.yaml to match the number of GPUs requested here.
#
# Usage
# -----
# From examples/openmm_cdk9_cyclinT1/:
#   sbatch scripts/hpc_we_apo.sl
#
# Adjust the #SBATCH directives below for your cluster.

#SBATCH --job-name=cdk9_we_apo
#SBATCH --output=logs/we_apo_%j.out
#SBATCH --error=logs/we_apo_%j.err
#SBATCH --time=48:00:00          # wall time limit (hh:mm:ss)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8        # one CPU per GPU for data I/O threads
#SBATCH --mem=32G
# Uncomment and edit for GPU nodes:
# #SBATCH --gres=gpu:4            # must match num_workers in config_apo.yaml
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

# ---------------------------------------------------------------------------
# Verify prerequisites
# ---------------------------------------------------------------------------
APO_BASIS="inputs/apo/basis_000.pdb"
if [ ! -f "$APO_BASIS" ]; then
    echo "ERROR: $APO_BASIS not found."
    echo "Run inputs/03_equilibrate.py first (sbatch scripts/hpc_equilibrate.sl)."
    exit 1
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p logs

cd "$(dirname "$0")/.."

echo "Starting CDK9 apo WE on $(hostname) at $(date)"
echo "SLURM job ID: ${SLURM_JOB_ID}"

python main.py --config config_apo.yaml

echo "CDK9 apo WE complete at $(date)"
