# CLAUDE.md — CDK9/CyclinT1 Phase 1 Weighted Ensemble

> **For AI coding assistants (Claude Code, Cursor, Copilot, etc.)**
> This file captures the full context of the CDK9 example — biology,
> architecture, decisions, current status, and how to continue development.
> Read this before reading any code.

---

## Biological context

**Target system:** CDK9 / CyclinT1 (P-TEFb transcription complex)
**Source structure:** PDB 4BCI — CDK9 (chain A) + CyclinT1 (chain B) +
  T3E inhibitor, resolved at 3.10 Å by X-ray crystallography.

**Scientific question:** Does CyclinT1 shift the CDK9 conformational
ensemble — specifically, does it lock the αC-helix in the active (αC-in)
state, or does CDK9 still sample αC-out conformations in the complex?

**Two conditions:**

| Condition | What's simulated | Expected behaviour |
|-----------|-----------------|-------------------|
| `apo` | CDK9 alone (chain A) | Free αC-helix mobility; both αC-in and αC-out accessible |
| `holo_cyclinT1` | CDK9 + CyclinT1 (chains A+B) | CyclinT1 expected to stabilise αC-in |

**T3E inhibitor:** removed from both conditions.
**Phospho-Thr186 (TPO):** retained — CDK9 is active (T-loop phosphorylated).
Requires AMBER ff14SB + PHOSAA10 force field parameters.

**CyclinT1 mutations Q77R/E96G/F241L:** all distal from CDK9 interface
(13.4 / 6.7 / 28.9 Å respectively).  Accepted as-is for Phase 1.
See `inputs/02_check_mutations.py` for the analysis.

---

## Progress coordinates (pcoord)

| Index | Observable | Range | Biology |
|-------|-----------|-------|---------|
| `pcoord[0]` | CDK9 Cα RMSD to reference (Å) | 0–6 Å (recycled above) | Global structural drift / unfolding guard |
| `pcoord[1]` | Glu66 Cδ – Lys48 Nζ distance (Å) | 3–20+ Å | αC-helix salt-bridge: ~3.5 Å = αC-in (active), >8 Å = αC-out (inactive) |

**Residue numbering** uses PDB 4BCI canonical CDK9 numbering (UniProt P50750).
Verify with `inputs/04_verify_pcoord_residues.py` before running.

---

## Architecture

This example subclasses four deepdrivewe extension points:

```
deepdrivewe base class              CDK9 subclass (this example)
─────────────────────               ────────────────────────────
ContactMapRMSDReporter     →    CDK9PcoordReporter    (simulate.py)
  adds salt-bridge distance; returns (n_frames, 2)

SimulationAgent            →    CDK9SimulationAgent   (simulate.py)
  injects CDK9PcoordReporter instead of ContactMapRMSDReporter

Recycler                   →    BoundaryRecycler      (recyclers.py)
  recycles walkers with pcoord[0] > rmsd_threshold (default 6 Å)

OrchestratorAgent          →    CDK9OrchestratorAgent (orchestrator.py)
  overrides evaluate_goals() with αC-in fraction logging

New class (no deepdrivewe base):
  Rectilinear2DBinner                                  (binners.py)
  2D uniform grid over (RMSD, salt-bridge distance)
```

Entry point: `main.py` — wires everything and calls `asyncio.run(run_academy_workflow(cfg))`.

---

## Phase roadmap

### Phase 1 — Geometric pcoord (THIS EXAMPLE, `feature/phase1-geometric-pcoord`)
**Status: code complete, awaiting HPC execution.**

- 2D pcoord: Cα RMSD + Glu66–Lys48 salt-bridge distance
- Uniform 2D RectilinearBinner (7 × 11 = 77 bins)
- BoundaryRecycler: RMSD > 6 Å → recycle to basis state
- HuberKim resampler (walkers merge/split per bin)
- Independent apo and holo_cyclinT1 runs
- Implicit solvent (GBn2) for portability

**Deliverable:** 2D probability flux maps over (RMSD, salt-bridge) for
apo vs holo; qualitative answer to whether CyclinT1 shifts the pcoord[1]
distribution.

### Phase 2 — CVAE latent pcoord (`feature/phase2-cvae-pcoord`, future)

- Replace geometric pcoord[1] with a CVAE latent coordinate
- Contact maps are already being collected (see `CDK9PcoordReporter`)
  so Phase 1 trajectories can seed Phase 2 training
- Voronoi binner in latent space
- Requires `deepdrivewe.academy_agents.training.TrainingAgent`

### Phase 3 — Comparative orchestration (`feature/phase3-comparative`, future)

- Coupled apo + holo ensembles in a single workflow
- Comparative probability flux analysis
- Reward-based allocation toward distinguishing regions

---

## Key design decisions

### Why 2D pcoord instead of 1D?

pcoord[0] (RMSD) is a *safety guard*, not a meaningful sampling coordinate.
It prevents walkers from drifting into globally unfolded states.  Without it,
a 1D salt-bridge binner would happily place walkers in bins where CDK9 has
unfolded — biologically meaningless and computationally wasteful.

The two dimensions are not independent: high RMSD tends to correlate with
large salt-bridge distance, but not always.  The 2D grid captures this.

### Why BoundaryRecycler (not LowRecycler)?

`LowRecycler` (deepdrivewe built-in) recycles walkers that fall *below* a
threshold — i.e., walkers that reach a target state.  CDK9 Phase 1 has no
target state: the goal is to *map* the landscape, not funnel toward a specific
conformation.  `BoundaryRecycler` recycles walkers that stray too *far*, which
prevents unfolding without biasing toward any conformation.

### Why implicit solvent (GBn2)?

Portability across compute environments: no periodic box, no pressure
coupling, no water equilibration.  Adequate for mapping which conformational
states exist and their rough relative probabilities.  Phase 2/3 should consider
explicit solvent for quantitative rate estimates.

### Why collect contact maps in Phase 1?

`CDK9PcoordReporter` inherits contact-map collection from
`ContactMapRMSDReporter` even though contact maps are not used in Phase 1.
This is deliberate: re-running simulations is expensive.  Phase 2 CVAE training
needs contact maps, and Phase 1 trajectories will seed Phase 2.

### Why not modify LowRecycler / RectilinearBinner directly?

We do not modify upstream deepdrivewe classes because:
1. This repo is a reference/dependency; upstream changes would be breaking.
2. Subclassing keeps the CDK9 code self-contained in this example directory.
3. Future contributors can replace CDK9 classes without touching deepdrivewe.

---

## How to continue development with an AI assistant

### Recommended startup prompt

When starting a new session on this codebase, give the AI this context:

> "I'm working on the CDK9/CyclinT1 Phase 1 weighted ensemble example in
> deepdrivewe (`examples/openmm_cdk9_cyclinT1/`).  The example is on the
> `feature/academy-agents` branch.  Read `examples/openmm_cdk9_cyclinT1/CLAUDE.md`
> and `CLAUDE.md` (repo root) first, then read `simulate.py` and `main.py`
> before making any changes.  The task is: [your task here]."

### Before making changes

1. Read `simulate.py` — most pcoord changes start here.
2. Read `main.py` — to understand how components wire together.
3. Check `inputs/04_verify_pcoord_residues.py` output — residue numbers
   in the structure may differ from defaults.

### Common tasks

**Add a 3rd pcoord dimension (e.g. DFG loop angle):**
1. Extend `CDK9PcoordReporter.report()` — append to `self._pcoord3`.
2. Update `get_rmsds()` to return `np.column_stack([rmsd, sb, dim3])`.
3. Replace `Rectilinear2DBinner` with a 3D binner or use `pcoord_idx=0`
   for the existing 1D binner.

**Change the recycling threshold:**
Edit `rmsd_boundary_ang` in `config_apo.yaml` / `config_holo_cyclinT1.yaml`.
No code change needed.

**Switch from implicit to explicit solvent:**
1. Change `solvent_type: explicit` in the YAML.
2. Add a `top_file` path (AMBER .prmtop).
3. Change `hardware_platform: CUDA` (explicit solvent needs GPU).

**Add Phase 2 CVAE analysis:**
1. Subclass `AnalysisPoolAgent` with a CVAE trainer.
2. Contact maps are already in `SimResult.data['contact_maps']`.
3. Replace `Rectilinear2DBinner` with a Voronoi binner over latent space.

---

## Known issues and open questions

- **Residue numbering:** Glu66/Lys48 numbers assume 4BCI PDB canonical
  numbering is preserved through PDBFixer.  Always verify with
  `inputs/04_verify_pcoord_residues.py` on a fresh structure.

- **holo system size:** ~4,700 atoms (vs ~2,600 for apo).  Budget memory
  accordingly when setting `num_workers`.

- **Implicit solvent drift:** GBn2 can cause unphysical loop extension in
  long segments.  The 6 Å RMSD boundary is conservative; reduce to 4 Å if
  you see frequent recycling.

- **Contact map shape mismatch:** apo (321 residues) and holo (573 residues)
  produce different contact map sizes.  Do not train a single CVAE on both
  conditions without alignment/padding.

- **PDB 4BCI duplicate chain IDs:** 4BCI stores protein and HETATM records
  for chain A and B, giving four chains total (A, B, A, B).  PDBFixer
  exposes this as duplicate IDs.  `inputs/01_download_and_clean.py` handles
  this with index-based chain removal — do not revert to name-based removal.

---

## How this example was developed

This example was built using **Claude Code** (claude-sonnet-4-6) as an agentic
coding assistant, following an iterative workflow:

1. **System context first** — the AI was given the biological question, PDB ID,
   and pcoord choices before any code was written.
2. **Read before write** — the AI read the deepdrivewe source classes
   (`ContactMapRMSDReporter`, `Recycler`, `OrchestratorAgent`, etc.) before
   subclassing them.
3. **Incremental commits** — each logical change was committed separately
   with descriptive messages, making the branch history self-explanatory.
4. **CLAUDE.md as living documentation** — this file was written as part of
   the implementation, not as an afterthought, so future AI sessions have
   full context from the start.

**Agentic development workflow for this codebase:**

```
Session start:
  → AI reads CLAUDE.md (this file) + repo CLAUDE.md
  → AI reads simulate.py + main.py
  → AI plans changes with the developer

Implementation:
  → Small, focused edits (not large rewrites)
  → Test with inputs/04_verify_pcoord_residues.py after structure changes
  → Commit each logical unit

Session end:
  → Update CLAUDE.md if new decisions were made
  → Commit with descriptive message so history is self-explanatory
```
