# Quickstart

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

python 01_download_and_clean.py          # ~1 min  — downloads PDB 4BCI, extracts apo/holo
python 02_check_mutations.py             # optional — reports CyclinT1 mutation distances
python 03_equilibrate.py                 # ~3–5 h CPU; submit via HPC script for GPU
python 04_verify_pcoord_residues.py      # verify Glu66/Lys48 residue numbering
```

On HPC, submit from `examples/openmm_cdk9_cyclinT1/`:

```bash
sbatch scripts/hpc_equilibrate.sl
```

## Step 2 — Run the weighted ensemble

```bash
# From examples/openmm_cdk9_cyclinT1/
export OPENMM_CPU_THREADS=1    # CPU only; remove this line for GPU

python main.py --config config_apo.yaml
python main.py --config config_holo_cyclinT1.yaml
```

On HPC:

```bash
sbatch scripts/hpc_we_apo.sl
sbatch scripts/hpc_we_holo.sl
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

## Memory budgeting

| Condition | Approx. atoms | Suggested `num_workers` |
|-----------|--------------|------------------------|
| apo | ~2,600 | 4–8 (GPU) |
| holo_cyclinT1 | ~4,700 | 2–4 (GPU) — larger memory footprint |

Adjust based on available GPU VRAM.  Each worker holds one OpenMM simulation
in memory simultaneously.
