# CDK9 / CyclinT1 Phase 1 — Geometric Pcoord Weighted Ensemble

An Academy-based weighted ensemble (WE) example for mapping the CDK9
αC-helix conformational landscape using a 2D geometric progress coordinate.

> **AI developers:** Read `CLAUDE.md` in this directory first — it has
> the full biological context, architecture overview, design decisions,
> and continuation guidance that CLAUDE.md is designed to provide.

---

## Scientific background

**System:** CDK9 (Ser/Thr kinase) / CyclinT1 (activating partner) — the
P-TEFb transcription elongation complex.  CDK9 is a therapeutic target
for cancer and viral infection (HIV-1 Tat).

**Question:** Does CyclinT1 shift CDK9's conformational ensemble, specifically
locking the αC-helix in the active (αC-in) state?

**Two conditions:**

| Condition | Simulation | Expected behaviour |
|-----------|-----------|-------------------|
| `apo` | CDK9 alone (PDB 4BCI chain A) | Free αC-helix mobility |
| `holo_cyclinT1` | CDK9 + CyclinT1 (chains A + B) | CyclinT1 expected to stabilise αC-in |

---

## Progress coordinates

| pcoord | Observable | Biology |
|--------|-----------|---------|
| `[0]` | CDK9 Cα RMSD to reference (Å) | Global structural drift; walkers recycled above 6 Å |
| `[1]` | Glu66 Cδ – Lys48 Nζ distance (Å) | αC-helix salt-bridge: ~3.5 Å = active, >8 Å = inactive |

---

## Example layout

```
examples/openmm_cdk9_cyclinT1/
├── CLAUDE.md                       AI assistant guidance (read first)
├── README.md                       This file
├── main.py                         Entry point
├── simulate.py                     CDK9PcoordReporter + CDK9SimulationAgent
├── recyclers.py                    BoundaryRecycler
├── binners.py                      Rectilinear2DBinner (2D grid)
├── orchestrator.py                 CDK9OrchestratorAgent
├── config_apo.yaml                 Apo condition configuration
├── config_holo_cyclinT1.yaml       Holo condition configuration
├── inputs/                         Structure preparation pipeline
│   ├── README.md
│   ├── 01_download_and_clean.py    Download 4BCI, prepare apo/holo PDBs
│   ├── 02_check_mutations.py       Inspect CyclinT1 mutations
│   ├── 03_equilibrate.py           Minimise + 2 ns NVT + save basis states
│   ├── 04_verify_pcoord_residues.py Verify Glu66/Lys48 numbering
│   └── .gitignore                  Ignores large/regenerable files
└── scripts/                        HPC job submission
    ├── README.md
    ├── hpc_equilibrate.sl          Slurm: structure equilibration
    ├── hpc_we_apo.sl               Slurm: CDK9 apo WE run
    └── hpc_we_holo.sl              Slurm: holo WE run
```

---

## Quickstart

### 1. Install deepdrivewe

```bash
# From the repo root
pip install -e .
```

### 2. Prepare structures

```bash
cd examples/openmm_cdk9_cyclinT1/inputs
python 01_download_and_clean.py          # ~1 min; downloads PDB 4BCI
python 02_check_mutations.py             # optional; checks CyclinT1 mutations
python 03_equilibrate.py                 # ~3-5 h on CPU; use HPC for GPU
python 04_verify_pcoord_residues.py      # verify residue numbering
```

For large clusters, submit `scripts/hpc_equilibrate.sl` instead.

### 3. Run the weighted ensemble

```bash
# From deepdrivewe repo root:
export OPENMM_CPU_THREADS=1  # CPU-only; remove for GPU

python examples/openmm_cdk9_cyclinT1/main.py \
    --config examples/openmm_cdk9_cyclinT1/config_apo.yaml

python examples/openmm_cdk9_cyclinT1/main.py \
    --config examples/openmm_cdk9_cyclinT1/config_holo_cyclinT1.yaml
```

For HPC: `sbatch examples/openmm_cdk9_cyclinT1/scripts/hpc_we_apo.sl`

### 4. Resume from checkpoint

Re-run the same command.  `EnsembleCheckpointer` automatically loads the
latest checkpoint.

---

## Configuration

Edit `config_apo.yaml` or `config_holo_cyclinT1.yaml`:

| Key | Default | Purpose |
|-----|---------|---------|
| `num_iterations` | 100 | WE iterations to run |
| `num_workers` | 4 | Parallel simulation workers (= GPUs) |
| `simulation_config.hardware_platform` | `CUDA` | `CUDA`, `OpenCL`, or `CPU` |
| `simulation_config.simulation_length_ns` | 0.05 | Segment length (50 ps) |
| `rmsd_boundary_ang` | 6.0 | RMSD recycling threshold (Å) |
| `sims_per_bin` | 4 | Target walkers per 2D bin cell |

---

## Expected outputs

```
runs/cdk9-apo/
├── params.yaml            Full resolved config
├── runtime.log            Structured application log
├── checkpoints/           Per-iteration ensemble state (HDF5)
└── simulations/
    └── iter_0001_sim_0000/
        ├── seg.pdb        Final frame → restart file for next iteration
        ├── seg.dcd        Trajectory
        ├── seg.log        Energy / temperature
        └── config.yaml    Simulation parameters
```

---

## Design notes

**Why 2D pcoord?**
pcoord[0] (RMSD) acts as a safety guard against global unfolding; it is not
meaningful for conformational analysis.  pcoord[1] (salt-bridge) is the
biologically relevant coordinate.  Combining them in a 2D grid keeps sampling
unbiased while preventing runaway trajectories.

**Why no target state?**
`BoundaryRecycler` replaces `LowRecycler`: we want to *map* the CDK9
landscape, not funnel toward a target.  Recycling happens only when walkers
stray too far (RMSD > 6 Å), not when they reach a specific conformation.

**Contact maps collected but unused in Phase 1:**
`CDK9PcoordReporter` still accumulates Cα contact maps (inherited from
`ContactMapRMSDReporter`).  These are stored in `SimResult.data['contact_maps']`
for Phase 2 CVAE training — avoiding the need to re-run simulations.

---

## Phased roadmap

| Phase | Branch | Status | Description |
|-------|--------|--------|-------------|
| 1 | `feature/phase1-geometric-pcoord` | **Code complete** | Geometric 2D pcoord, independent conditions |
| 2 | `feature/phase2-cvae-pcoord` | Planned | CVAE latent pcoord, Voronoi binner |
| 3 | `feature/phase3-comparative` | Planned | Coupled apo+holo comparison, adaptive guidance |

---

## Comparison with NTL9 Academy example

| Feature | NTL9 (`openmm_ntl9_hk_academy`) | CDK9 (`openmm_cdk9_cyclinT1`) |
|---------|--------------------------------|-------------------------------|
| pcoord dims | 1 (RMSD) | 2 (RMSD + salt-bridge) |
| Binner | 1D `RectilinearBinner` | 2D `Rectilinear2DBinner` |
| Recycler | `LowRecycler` (target-directed) | `BoundaryRecycler` (no target) |
| Reporter | `ContactMapRMSDReporter` | `CDK9PcoordReporter` (subclass) |
| Agent | `SimulationAgent` | `CDK9SimulationAgent` (subclass) |
| Goal loop | Placeholder | αC-in fraction logging |
