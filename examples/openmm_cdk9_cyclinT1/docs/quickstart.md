# Quickstart

## Compute expectations

The example has two distinct runtime profiles:

| Stage | Local CPU | GPU (CUDA) | Notes |
|-------|-----------|------------|-------|
| `01_download_and_clean.py` | ~1 min | — | Downloads PDB 4BCI |
| `02_check_mutations.py` | ~1 min | — | Optional; informational only |
| `03_equilibrate.py` | **3–5 h** | ~20–40 min | Both conditions; long-running |
| `04_verify_pcoord_residues.py` | ~1 min | — | Must pass before WE run |
| WE run (100 iter, 4 workers) | Not practical | Hours–days | Designed for multi-GPU |

**`03_equilibrate.py` is the critical long-running step.** Plan accordingly —
run it on a GPU workstation or overnight on CPU.  The WE workflow is designed
for multi-GPU execution; CPU-only runs are useful for smoke-testing the
pipeline (2 iterations, 1 ps segments) but not for production sampling.

---

## Prerequisites

```bash
# From deepdrivewe repo root
pip install -e .

# OpenMM + PDBFixer via conda-forge if not already installed:
# conda install -c conda-forge openmm pdbfixer
```

## Step 1 — Prepare structures

```bash
cd examples/openmm_cdk9_cyclinT1/inputs

python 01_download_and_clean.py     # ~1 min — downloads PDB 4BCI, extracts apo/holo PDBs
python 02_check_mutations.py        # optional — reports CyclinT1 mutation distances
python 03_equilibrate.py            # 3–5 h CPU or ~30 min GPU — plan accordingly
```

Before proceeding, verify the equilibrated structures are usable:

```bash
python 04_verify_pcoord_residues.py
```

This script checks that Glu66 and Lys48 are present with the expected atom
names in both prepared structures.  **Do not proceed to Step 2 if this fails**
— a residue numbering shift from PDBFixer will silently produce wrong pcoord
values throughout the WE run.

Expected output of a passing check:
```
apo:  Glu66 CD found at index N  |  Lys48 NZ found at index M
holo: Glu66 CD found at index N  |  Lys48 NZ found at index M
All checks passed.
```

## Step 2 — Run the weighted ensemble

```bash
# From examples/openmm_cdk9_cyclinT1/
export OPENMM_CPU_THREADS=1    # CPU only; remove this line for GPU

python main.py --config config_apo.yaml
python main.py --config config_holo_cyclinT1.yaml
```

## Step 3 — Resume from checkpoint

Re-run the same command — `EnsembleCheckpointer` loads the latest checkpoint
automatically.

---

## Configuration reference

Edit `config_apo.yaml` or `config_holo_cyclinT1.yaml`:

| Key | Default | Purpose |
|-----|---------|---------|
| `num_iterations` | 100 | WE iterations to run |
| `num_workers` | 4 | Parallel simulation workers (= GPUs for GPU runs) |
| `simulation_config.hardware_platform` | `CUDA` | `CUDA`, `OpenCL`, or `CPU` |
| `simulation_config.simulation_length_ns` | 0.05 | Segment length (50 ps) |
| `rmsd_boundary_ang` | 6.0 | RMSD recycling threshold (Å) — reduce to 4.0 if frequent recycling |
| `sims_per_bin` | 4 | Target walkers per 2D bin cell |

### CPU smoke test

Temporarily edit `config_apo.yaml`:

```yaml
num_iterations: 2
num_workers: 1
simulation_config:
  hardware_platform: CPU
  simulation_length_ns: 0.001
```

```bash
export OPENMM_CPU_THREADS=1
python main.py --config config_apo.yaml
```

---

## Expected outputs

```
runs/cdk9-apo/
├── params.yaml              Full resolved config
├── runtime.log              Application log
├── checkpoints/             Per-iteration ensemble state (HDF5)
└── simulations/
    └── iter_0001_sim_0000/
        ├── seg.pdb          Final frame — restart file for next iteration
        ├── seg.dcd          Trajectory
        ├── seg.log          Energy / temperature
        └── config.yaml      Simulation parameters
```

---

## Running on an HPC cluster (optional)

The example runs with `python main.py` on any machine with OpenMM installed.
The `scripts/` directory contains Slurm job scripts as a **starting-point
template** for HPC submission — they are not required and will not work
without cluster-specific edits.

Before submitting any script you must:

1. Add your environment activation (e.g. `conda activate myenv` or
   `source .venv/bin/activate`) in the `# Environment setup` section
2. Uncomment and set `#SBATCH --gres=gpu:N` and `--partition` to match
   your cluster's GPU resource syntax
3. Set `num_workers` in the relevant YAML config to match the GPU count

See `scripts/README.md` for the full checklist.

---

## Memory budgeting

| Condition | Approx. atoms | Suggested `num_workers` |
|-----------|--------------|------------------------|
| apo | ~2,600 | 4–8 (GPU) |
| holo_cyclinT1 | ~4,700 | 2–4 (GPU) — larger memory footprint |

Adjust based on available GPU VRAM.  Each worker holds one OpenMM simulation
in memory simultaneously.
