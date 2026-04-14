# CLAUDE.md — CDK9 / CyclinT1 Phase 1 Weighted Ensemble

> **For AI coding assistants (Claude Code, Cursor, Copilot, etc.)**
>
> This is the single reference document for this example.  Read it before
> touching any code.  It covers the biology, the deepdrivewe architecture,
> every design decision made, the phased roadmap, and explicit instructions
> for how to continue AI-assisted development on this codebase.

---

## 1. What this example is

A complete, self-contained weighted ensemble (WE) simulation example for
mapping the CDK9 / CyclinT1 conformational landscape.  It lives at:

```
deepdrivewe/examples/openmm_cdk9_cyclinT1/
```

It is built on the **deepdrivewe + Academy agents** framework and follows the
same pattern as `examples/openmm_ntl9_hk_academy/`, extended with a 2D
geometric progress coordinate tailored to CDK9 biology.

**Key documents for this session:**
- `dev/intent.md` — scientific goals and biological motivation
- `dev/spec.md` — technical design decisions and decision log
- `dev/phases/phase1.md` — Phase 1 plan and status (code complete)
- `docs/quickstart.md` — how to run the example

---

## 2. Biological context

**Target:** CDK9 (Ser/Thr kinase) / CyclinT1 (activating partner) — the
P-TEFb transcription elongation complex.  CDK9 is a cancer and antiviral drug
target.

**Source structure:** PDB 4BCI — CDK9 (chain A) + CyclinT1 (chain B) +
T3E inhibitor, 3.10 Å resolution.

**Scientific question:** Does CyclinT1 shift CDK9's conformational ensemble
— specifically, does it lock the αC-helix in the active (αC-in) state?

| Condition | What runs | Expected |
|-----------|----------|---------|
| `apo` | CDK9 alone (chain A) | Free αC-helix mobility: both αC-in and αC-out sampled |
| `holo_cyclinT1` | CDK9 + CyclinT1 (chains A+B) | CyclinT1 stabilises αC-in |

**Key structural notes:**
- T3E inhibitor removed from both conditions (we study the unliganded landscape).
- Phospho-Thr186 (TPO) retained — active CDK9 is Thr186-phosphorylated.
  Requires AMBER ff14SB + PHOSAA10 force field parameters.
- CyclinT1 mutations Q77R/E96G/F241L: all distal (13.4/6.7/28.9 Å from CDK9
  interface); accepted as-is for Phase 1 (see `inputs/02_check_mutations.py`).
- PDB 4BCI has duplicate chain IDs (A, B, A, B) — second pair is HETATM.
  Chain removal uses index-based logic in `inputs/01_download_and_clean.py`.

---

## 3. Progress coordinates

| Index | Observable | Range | Biology |
|-------|-----------|-------|---------|
| `pcoord[0]` | CDK9 Cα RMSD to reference (Å) | 0–6 Å (recycled above) | Global structural drift / unfolding guard |
| `pcoord[1]` | Glu66 Cδ – Lys48 Nζ distance (Å) | 3–20+ Å | αC-helix salt-bridge: ~3.5 Å = αC-in (active), >8 Å = αC-out |

Residue numbering follows PDB 4BCI / UniProt P50750.
Verify with `inputs/04_verify_pcoord_residues.py` before running.

---

## 4. deepdrivewe architecture (what you need to know)

deepdrivewe provides the WE framework.  This example subclasses four extension
points; everything else is inherited unchanged.

### 4.1 Package layout

```
deepdrivewe/deepdrivewe/          ← Python package (importable)
├── api.py                        SimMetadata, BasisStates, TargetState
├── academy_agents/
│   ├── orchestrator.py           OrchestratorAgent — workflow loop
│   ├── ensemble.py               EnsembleManagerAgent — binner/resampler/recycler
│   ├── simulation.py             SimulationAgent — runs one OpenMM segment
│   └── config.py                 SimulationPoolConfig, AcademyWorkflowConfig
├── simulation/
│   └── openmm.py                 OpenMMConfig, OpenMMSimulation,
│                                 ContactMapRMSDReporter
├── binners/
│   └── rectilinear.py            RectilinearBinner (1D only)
├── recyclers/
│   └── base.py                   Recycler ABC + recycle_simulations()
└── resamplers/                   HuberKimResampler
```

### 4.2 CDK9 subclasses (this example)

```
ContactMapRMSDReporter → CDK9PcoordReporter    (simulate.py)
SimulationAgent        → CDK9SimulationAgent   (simulate.py)
Recycler               → BoundaryRecycler      (recyclers.py)
Binner (new)           → Rectilinear2DBinner   (binners.py)
OrchestratorAgent      → CDK9OrchestratorAgent (orchestrator.py)
```

### 4.3 Key data flow

```
CDK9SimulationAgent.run_simulation()
  → CDK9PcoordReporter.report() on every frame
      → parent accumulates contact maps (for Phase 2 CVAE)
      → child accumulates salt-bridge distance
  → CDK9PcoordReporter.get_rmsds() → (n_frames, 2) array
  → stored in SimMetadata.pcoord (list[list[float]])

EnsembleManagerAgent.apply_resampling()
  → Rectilinear2DBinner.assign_bins(pcoords)
      pcoords[:, 0] = last-frame RMSD per walker
      pcoords[:, 1] = last-frame salt-bridge dist per walker
  → BoundaryRecycler.recycle(pcoords) → indices where RMSD > threshold
  → HuberKimResampler merges/splits walkers within each bin
```

### 4.4 pcoord conventions

- `SimMetadata.pcoord` shape: `(n_frames, n_dims)` — stored as list of lists.
- `ContactMapRMSDReporter.get_rmsds()` returns `(n_frames, 1)` by default.
- `CDK9PcoordReporter.get_rmsds()` returns `(n_frames, 2)`.
- Recycler and binner receive `pcoords[-1]` — last frame of each walker.

### 4.5 Academy agent framework

Academy uses async Python with four key primitives:

```python
from academy.agent import action, loop
from academy.handle import Handle
from academy.manager import Manager
from academy.exchange import LocalExchangeFactory

# @action  — callable remotely via Handle.method_name(args)
# @loop    — runs in background until shutdown event is set
# Manager  — launches agents, creates handles
# Handle   — typed proxy to a remote agent
```

Agents communicate via handles only (no direct calls).  The `@action`
decorator must be preserved on overridden methods or the framework won't
register them.

---

## 5. Design decisions and rationale

### Why 2D pcoord instead of 1D?
pcoord[0] (RMSD) is a *safety guard*, not a sampling coordinate.  Without it,
a 1D binner over salt-bridge distance alone would keep walkers in biologically
meaningless unfolded states.  The two dimensions are complementary, not
redundant.

### Why BoundaryRecycler, not LowRecycler?
`LowRecycler` (deepdrivewe built-in) recycles walkers that fall *below* a
threshold — i.e., walkers reaching a target state.  CDK9 Phase 1 has no target
state; the goal is to *map* the landscape.  `BoundaryRecycler` recycles walkers
that stray too *far* (RMSD > 6 Å), preventing unfolding without biasing toward
any specific conformation.

### Why implicit solvent (GBn2)?
Portability: no periodic box, no pressure coupling, no water equilibration.
Adequate for mapping which states exist and their rough probabilities.  Phase
2/3 should consider explicit solvent for quantitative rate estimates.

### Why collect contact maps in Phase 1 if they're unused?
`CDK9PcoordReporter` inherits contact-map collection from the parent.  Phase 2
CVAE training needs contact maps, and Phase 1 trajectories will seed Phase 2.
Re-running simulations is expensive; collect now, use later.

### Why subclass instead of modifying deepdrivewe base classes?
deepdrivewe is a reference library; upstream modifications are breaking changes
for all other examples.  CDK9 code is self-contained in this directory and can
be replaced/extended without touching deepdrivewe.

---

## 6. Example layout

```
examples/openmm_cdk9_cyclinT1/
├── CLAUDE.md                        ← This file (AI assistant guide)
├── README.md                        ← Broad overview and navigation
├── docs/
│   └── quickstart.md                ← How to run: prereqs, steps, config, outputs
├── dev/
│   ├── intent.md                    ← Scientific goals and biological motivation
│   ├── spec.md                      ← Technical design decisions and decision log
│   └── phases/
│       ├── phase1.md                ← Phase 1 plan and deliverables (complete)
│       ├── phase2.md                ← Phase 2 plan: CVAE latent pcoord (planned)
│       └── phase3.md                ← Phase 3 plan: comparative orchestration (planned)
├── main.py                          ← Entry point; wires Academy agents
├── simulate.py                      ← CDK9PcoordReporter + CDK9SimulationAgent
├── recyclers.py                     ← BoundaryRecycler
├── binners.py                       ← Rectilinear2DBinner (2D uniform grid)
├── orchestrator.py                  ← CDK9OrchestratorAgent
├── config_apo.yaml                  ← Full config for apo condition
├── config_holo_cyclinT1.yaml        ← Full config for holo condition
├── inputs/
│   ├── README.md                    ← Pipeline docs + structural decisions table
│   ├── 01_download_and_clean.py     ← Download 4BCI, prepare apo/holo PDBs
│   ├── 02_check_mutations.py        ← CyclinT1 mutation proximity analysis
│   ├── 03_equilibrate.py            ← Minimise + 2 ns NVT + save basis states
│   ├── 04_verify_pcoord_residues.py ← Verify Glu66/Lys48 numbering
│   └── .gitignore
└── scripts/
    ├── README.md                    ← HPC submission guide
    ├── hpc_equilibrate.sl           ← Slurm: structure equilibration
    ├── hpc_we_apo.sl                ← Slurm: CDK9 apo WE run
    └── hpc_we_holo.sl               ← Slurm: holo WE run
```

---

## 7. Running the example

### Prerequisites

```bash
# From deepdrivewe repo root
pip install -e .                       # install deepdrivewe + dependencies
# openmm via conda-forge if not present:
# conda install -c conda-forge openmm pdbfixer
```

### Step 1 — Prepare structures

```bash
cd examples/openmm_cdk9_cyclinT1/inputs
python 01_download_and_clean.py          # ~1 min
python 02_check_mutations.py             # informational
python 03_equilibrate.py                 # ~3-5 h CPU; use HPC for GPU
python 04_verify_pcoord_residues.py      # verify Glu66/Lys48
```

Or on HPC: `sbatch scripts/hpc_equilibrate.sl` (from within this directory).

### Step 2 — Run the workflow

```bash
# From this directory (examples/openmm_cdk9_cyclinT1/):
export OPENMM_CPU_THREADS=1    # CPU only; remove for GPU

python main.py --config config_apo.yaml
python main.py --config config_holo_cyclinT1.yaml
```

On HPC: `sbatch scripts/hpc_we_apo.sl` / `sbatch scripts/hpc_we_holo.sl`.

### Step 3 — Resume from checkpoint

Re-run the same command.  `EnsembleCheckpointer` loads the latest checkpoint
automatically.

---

## 8. Common tasks for AI assistants

### Changing the recycling threshold
Edit `rmsd_boundary_ang` in the YAML config.  No code change needed.

### Adding a 3rd pcoord dimension
1. Extend `CDK9PcoordReporter.report()` in `simulate.py` — append to a
   new `self._pcoord3` list.
2. Update `get_rmsds()` to return `np.column_stack([rmsd, sb, dim3])`.
3. Replace `Rectilinear2DBinner` with a 3D variant or use the 1D
   `RectilinearBinner` on a single dimension.

### Switching to explicit solvent
1. Set `solvent_type: explicit` in the YAML.
2. Add a `top_file` path (AMBER .prmtop).
3. Set `hardware_platform: CUDA`.

### Adding Phase 2 CVAE analysis
1. Contact maps are already in `SimResult.data['contact_maps']`.
2. Subclass `AnalysisPoolAgent` with a CVAE trainer.
3. Replace `Rectilinear2DBinner` with a Voronoi binner over latent space.

### Verifying the example runs end-to-end (CPU smoke test)

```bash
# Edit config_apo.yaml temporarily:
#   num_iterations: 2
#   num_workers: 1
#   simulation_config.hardware_platform: CPU
#   simulation_config.simulation_length_ns: 0.001
export OPENMM_CPU_THREADS=1
python main.py --config config_apo.yaml
```

---

## 9. Known issues and open questions

- **Residue numbering:** Glu66/Lys48 assume 4BCI canonical CDK9 numbering
  is preserved through PDBFixer.  Always verify with
  `inputs/04_verify_pcoord_residues.py` on a fresh structure.

- **holo system size:** ~4,700 atoms vs ~2,600 for apo.  Budget memory
  accordingly when setting `num_workers`.

- **GBn2 loop drift:** Implicit solvent can cause unphysical loop extension
  in long segments.  Reduce `rmsd_boundary_ang` to 4 Å if frequent recycling
  is observed.

- **Contact map shape mismatch:** apo and holo produce different contact map
  sizes (321 vs 573 residues).  Do not train a single CVAE on both conditions
  without alignment/padding.

- **`@action` on overridden methods:** If you override any `@action`-decorated
  method from a deepdrivewe base class, keep the `@action` decorator — the
  Academy framework won't register it otherwise.

---

## 10. How this example was developed (agentic workflow)

This example was built using **Claude Code** (claude-sonnet-4-6) as an
agentic coding assistant.  The development process is captured here so
future AI-assisted sessions can start with full context rather than
re-deriving it.

### What was done

1. **System context first** — biological question, PDB ID, pcoord design, and
   structural decisions were established before any code was written.
2. **Read before write** — the AI read deepdrivewe source classes
   (`ContactMapRMSDReporter`, `Recycler`, `SimulationAgent`, etc.) before
   subclassing them.
3. **Self-contained example** — all CDK9-specific code lives in this directory;
   nothing was added to deepdrivewe's core package.
4. **CLAUDE.md written as part of the implementation** — not as an afterthought,
   so future AI sessions have full context from the start.

### Recommended startup prompt for future sessions

> "I'm working on the CDK9/CyclinT1 WE example at
> `examples/openmm_cdk9_cyclinT1/` in the deepdrivewe repo
> (`feature/academy-agents` branch).  Read `CLAUDE.md` in that directory
> first, then `dev/intent.md`, `dev/spec.md`, and the relevant
> `dev/phases/` file before making any changes.
> The task is: [your task here]."

### Session workflow

```
Session start:
  → AI reads CLAUDE.md (this file)
  → AI reads simulate.py and main.py
  → AI plans changes with the developer before writing code

Implementation:
  → Small, focused edits — not large rewrites
  → Test with 04_verify_pcoord_residues.py after any structural changes
  → Commit each logical unit with a descriptive message

Session end:
  → Update CLAUDE.md if new decisions were made
  → Update ROADMAP.md if phase status changed
  → Commit so history is self-explanatory for the next session
```
