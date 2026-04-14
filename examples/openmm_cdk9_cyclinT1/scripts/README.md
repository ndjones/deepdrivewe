# scripts/ — HPC submission scripts

Generic Slurm scripts for running the CDK9 Phase 1 pipeline on an HPC cluster.

## Workflow order

```
1. hpc_equilibrate.sl   Minimise + equilibrate, save basis states
2. hpc_we_apo.sl        CDK9 apo weighted ensemble
   hpc_we_holo.sl       CDK9/CyclinT1 holo weighted ensemble  (can run in parallel with apo)
```

## Before submitting

1. **Activate your Python environment** — edit the "Environment setup" section
   in each script for your cluster (conda, venv, or module-based).

2. **Configure GPU/partition directives** — the `#SBATCH --gres=gpu:N` and
   `--partition` lines are commented out.  Uncomment and set to match your
   cluster's GPU resource syntax.

3. **Match `num_workers` to GPU count** — `num_workers` in
   `cdk9_cyclinT1/config_apo.yaml` (and `config_holo_cyclinT1.yaml`)
   must equal the number of GPUs requested in the Slurm script.

4. **Set wall time** — default is 48 h for WE jobs.  Adjust to your
   cluster's partition limits.

## Submitting

```bash
# From the project root:
sbatch scripts/hpc_equilibrate.sl

# After basis states are ready:
sbatch scripts/hpc_we_apo.sl
sbatch scripts/hpc_we_holo.sl
```

## Logs

Job stdout/stderr go to `logs/` (created automatically).
Application logs go to `runs/<condition>/runtime.log`.

## CPU-only runs

If no GPUs are available, set `hardware_platform: CPU` in the YAML config
and remove `--gres=gpu` from the Slurm script.  Expect roughly 10–20×
slower simulation throughput compared to CUDA.
