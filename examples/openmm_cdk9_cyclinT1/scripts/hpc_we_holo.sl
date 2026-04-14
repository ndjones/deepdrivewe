#!/bin/bash
# Slurm job script — CDK9/CyclinT1 holo weighted ensemble (Phase 1)
#
# Runs the CDK9+CyclinT1 holo weighted ensemble using Academy agents.
# Each CDK9SimulationAgent worker occupies one GPU.  Set num_workers in
# config_holo_cyclinT1.yaml to match the GPUs requested here.
#
# Usage
# -----
# From examples/openmm_cdk9_cyclinT1/:
#   sbatch scripts/hpc_we_holo.sl
#
# Adjust the #SBATCH directives below for your cluster.

#SBATCH --job-name=cdk9_we_holo
#SBATCH --output=logs/we_holo_%j.out
#SBATCH --error=logs/we_holo_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G                # holo system is larger (~4700 atoms vs ~2600)
# Uncomment and edit for GPU nodes:
# #SBATCH --gres=gpu:4            # must match num_workers in config_holo_cyclinT1.yaml
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
HOLO_BASIS="inputs/holo_cyclinT1/basis_000.pdb"
if [ ! -f "$HOLO_BASIS" ]; then
    echo "ERROR: $HOLO_BASIS not found."
    echo "Run inputs/03_equilibrate.py first (sbatch scripts/hpc_equilibrate.sl)."
    exit 1
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
mkdir -p logs

cd "$(dirname "$0")/.."

echo "Starting CDK9/CyclinT1 holo WE on $(hostname) at $(date)"
echo "SLURM job ID: ${SLURM_JOB_ID}"

python main.py --config config_holo_cyclinT1.yaml

echo "CDK9/CyclinT1 holo WE complete at $(date)"
