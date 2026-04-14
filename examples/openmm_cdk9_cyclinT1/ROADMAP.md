# ROADMAP — CDK9 / CyclinT1 Weighted Ensemble Example

> **For AI coding assistants:** Read `CLAUDE.md` first.  This roadmap covers
> the scientific intent, phased plan, decision log, and open questions.  It is
> the source of truth for what has been done and what comes next.

---

## Scientific intent

This example adapts the deepdrivewe + Academy weighted ensemble (WE) framework
from its reference use case (protein folding) to a new goal: **mapping the
conformational landscape of CDK9 and understanding how its obligate activating
partner CyclinT1 shifts that landscape**, with the aim of explaining drug
selectivity mechanisms.

### Biological motivation

CDK9 is a serine/threonine kinase that, when complexed with CyclinT1, forms
P-TEFb (Positive Transcription Elongation Factor b).  Drug potency against CDK9
is modulated by the presence of CyclinT1 in ways that are not fully understood
mechanistically.  The hypothesis is that CyclinT1 alters which conformational
states of CDK9 are accessible — particularly the position of the αC-helix and
the DFG loop — and that this shift controls which inhibitor binding modes are
favoured, thereby driving selectivity differences.

To test this we need to:

1. Sample the conformational landscape of **CDK9 alone** (apo — without CyclinT1)
2. Sample the landscape of **CDK9 bound to CyclinT1** (holo)
3. Quantitatively compare the two landscapes to identify partner-driven shifts
4. Correlate those shifts with known selectivity and potency differences across
   inhibitor series

Weighted ensemble simulation is well-suited here because it provides
statistically rigorous estimates of the probability weight assigned to each
region of conformational space — not just which states exist, but how likely
they are.

---

## System: PDB 4BCI

**Structure:** CDK9 (chain A, 331 residues) + CyclinT1 (chain B, 260 residues)
\+ inhibitor T3E (chain C) from *Homo sapiens*, resolved at 3.10 Å by X-ray
diffraction.

**For simulation:**
- **Remove:** T3E ligand — studying the unliganded landscape
- **Keep (holo):** CDK9 (chain A) + CyclinT1 (chain B)
- **CDK9 only (apo):** chain A alone

| Feature | Detail | Action |
|---------|--------|--------|
| TPO (phosphoThr186) | Activation loop, phosphorylated in crystal | Retain — reflects active CDK9; requires AMBER ff14SB + PHOSAA10 |
| CyclinT1 mutations | Q77R/E96G/F241L in chain B (crystallisation construct) | All distal (13–29 Å from interface); accepted as-is for Phase 1 |
| Resolution 3.10 Å | Moderate; some loops poorly resolved | PDBFixer used to add missing atoms/hydrogens |
| Duplicate chain IDs | PDB has A/B/A/B — second pair is HETATM | Index-based removal in `inputs/01_download_and_clean.py` |

---

## Phased plan

---

### Phase 1 — Geometric pcoord, independent conditions

**Branch:** `feature/phase1-geometric-pcoord`
**Status:** Code complete

**Intent:**
Establish the minimum viable landscape-sampling workflow.  Replace the
folding-oriented single-scalar RMSD pcoord with a 2D geometric pcoord grounded
in CDK9 biology, remove the target-directed recycler, and run apo vs. holo as
two independent experiments compared post-hoc.  Validates the pipeline on the
real system before introducing ML-based pcoords.

No environment-specific compute optimisations in this phase.  The workflow is
portable across desktop CPU, HPC cluster, and cloud GPU; `num_workers` and
simulation length are the only knobs that need tuning per environment.

#### 1.1 System preparation

Two conditions, each requiring equilibrated basis states:

**Condition A — apo CDK9**
- Extract chain A from 4BCI; remove T3E + CyclinT1
- Retain phospho-Thr186 (TPO) — AMBER ff14SB + PHOSAA10
- GBn2 implicit solvent; minimise + 2 ns NVT equilibration
- Save 10 snapshots as basis states (`inputs/apo/*.pdb`)

**Condition B — holo CDK9 + CyclinT1**
- Extract chains A + B; remove T3E only
- Same phospho-Thr186 treatment
- GBn2 implicit solvent; minimise + 2 ns NVT equilibration
- Save 10 snapshots as basis states (`inputs/holo_cyclinT1/*.pdb`)

**Solvent choice:** Implicit solvent (GBn2) for Phase 1 — reduces compute cost
and is acceptable for initial landscape exploration.  Phase 2/3 should move to
explicit solvent for publication-grade rates.

#### 1.2 Progress coordinate

| Index | Observable | Biology |
|-------|-----------|---------|
| `pcoord[0]` | CDK9 Cα RMSD to reference (Å) | Global drift guard; recycle above 6 Å |
| `pcoord[1]` | Glu66 Cδ – Lys48 Nζ distance (Å) | αC-helix salt bridge; ~3.5 Å = αC-in (active), >8 Å = αC-out |

Glu66–Lys48 is the canonical kinase activation metric.  CyclinT1 stabilises
αC-in; apo CDK9 is expected to sample αC-out more frequently.

**Alternatives considered and set aside:**

| Option | Metric | Reason deferred |
|--------|--------|----------------|
| B | Phe168 – gatekeeper (Phe103) distance | DFG-in/out; less directly coupled to CyclinT1 |
| C | αC-helix RMSD (residues ~58–75) | Coarser; harder to interpret than the salt bridge |
| D | CDK9–CyclinT1 interface RMSD | Undefined for apo condition |

Options B/C/D remain available for Phase 2 if the salt bridge proves
insufficient.

#### 1.3 Binning

2D uniform grid (7 RMSD bins × 11 salt-bridge bins = 77 cells) via the new
`Rectilinear2DBinner`.  Flat index = `row * n_cols + col`.

#### 1.4 Recycling

`BoundaryRecycler` replaces `LowRecycler`.  Recycles walkers where `pcoord[0]`
(RMSD) exceeds `rmsd_boundary_ang` (default 6 Å), preventing unfolded states
from polluting the ensemble.  No target state — the goal is to *map* the
landscape, not funnel toward a specific conformation.

#### 1.5 Comparative analysis (post-hoc)

Run two independent workflows and compare:
- 2D bin occupation heatmaps (apo vs. holo, plus difference map)
- Marginal weight distribution over pcoord[1] — the key quantitative result
- Contact map differences: per-residue-pair contact probability apo vs. holo
- CVAE latent space overlay (Phase 2 warm-start)

Notebook: `analysis/compare_landscapes.ipynb` (planned; not included in Phase 1)

#### 1.6 `evaluate_goals` — coverage tracking

`CDK9OrchestratorAgent` overrides `evaluate_goals` to log per-iteration:
- Fraction of walkers with pcoord[1] < 5 Å (αC-in fraction)
- Per-iteration log message via `logger.info`

#### Phase 1 deliverables

- [x] `inputs/01_download_and_clean.py` — download 4BCI, prepare apo/holo PDBs
- [x] `inputs/02_check_mutations.py` — CyclinT1 mutation proximity analysis
- [x] `inputs/03_equilibrate.py` — minimise + 2 ns NVT + save basis states
- [x] `inputs/04_verify_pcoord_residues.py` — verify Glu66/Lys48 numbering
- [x] `BoundaryRecycler` — `recyclers.py`
- [x] `CDK9PcoordReporter` + `CDK9SimulationAgent` — `simulate.py`
- [x] `Rectilinear2DBinner` — `binners.py`
- [x] `CDK9OrchestratorAgent` — `orchestrator.py`
- [x] `config_apo.yaml` and `config_holo_cyclinT1.yaml`
- [x] `main.py` entry point
- [x] HPC Slurm scripts — `scripts/`
- [x] `CLAUDE.md` + `ROADMAP.md` — AI assistant context

---

### Phase 2 — CVAE latent pcoord, adaptive binning

**Branch:** `feature/phase2-cvae-pcoord`
**Status:** Planned — depends on Phase 1 trajectory data

**Intent:**
Replace the hand-crafted 2D geometric pcoord with CVAE latent coordinates
learned from contact maps collected in Phase 1.  Removes the need to
pre-specify which geometric features matter.

**Key changes:**
- Warm-start CVAE on Phase 1 contact maps before beginning Phase 2 WE run
- Wire `TrainingAgent` + `InferenceAgent` into the CDK9 workflow
- Replace `Rectilinear2DBinner` with a `VoronoiBinner` over latent space

**Contact map definition:**
- Intra-CDK9 Cα contacts only (keeps feature space consistent between apo/holo)
- Optionally add CDK9–CyclinT1 interface contacts as a second channel
- `cutoff_angstrom`: 8.0 Å (existing default)

#### Phase 2 deliverables

- [ ] CVAE warm-start training script (Phase 1 contact maps)
- [ ] `TrainingAgent` + `InferenceAgent` wired into CDK9 example
- [ ] `VoronoiBinner` in `deepdrivewe/binners/`
- [ ] Updated CDK9 configs for CVAE-pcoord runs
- [ ] Updated analysis notebook: latent space comparison apo vs. holo

---

### Phase 3 — Coupled comparative ensembles

**Branch:** `feature/phase3-comparative-orchestration`
**Status:** Planned — depends on Phase 2 CVAE latent pcoord

**Intent:**
Extend the orchestrator to manage apo and holo ensembles within a single
workflow, enabling online comparison.  Adaptively focuses compute on
conformational regions where the two CDK9 landscapes diverge most.

**Key changes:**
- `ComparativeOrchestratorAgent`: manages two `EnsembleManagerAgent` pairs;
  each iteration advances both, computes Jensen-Shannon divergence between bin
  weight distributions, allocates more walkers to diverging regions
- Divergence-driven resampler: upweights high-divergence bins
- Shared CVAE trained on pooled apo + holo trajectories with balanced batches

#### Phase 3 deliverables

- [ ] `ComparativeOrchestratorAgent`
- [ ] Divergence metric computation between ensembles
- [ ] Divergence-driven resampler
- [ ] Shared CVAE training on pooled data
- [ ] Updated analysis notebook with online divergence tracking

---

## Decision log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-04-14 | Target system: CDK9/CyclinT1, PDB 4BCI | Defined by biological question; T3E removed, phospho-Thr186 retained |
| 2026-04-14 | Binding partner: CyclinT1 (chain B in 4BCI) | Obligate activating subunit; directly implicated in inhibitor selectivity |
| 2026-04-14 | pcoord[1]: Glu66–Lys48 salt bridge | Canonical CDK activation metric; directly reports CyclinT1 effect; linked to inhibitor selectivity |
| 2026-04-14 | Start with Phase 1 geometric pcoord | No existing MD data; avoids CVAE cold-start; validates pipeline first |
| 2026-04-14 | Apo vs. holo as independent runs in Phase 1 | Simpler than coupled orchestration; sufficient for initial landscape comparison |
| 2026-04-14 | Implicit solvent (GBn2) for Phase 1 | Reduces compute; acceptable for landscape topology; explicit solvent flagged for Phase 2+ |
| 2026-04-14 | BoundaryRecycler instead of LowRecycler | No target state; recycle on RMSD upper bound only to prevent unfolding |
| 2026-04-14 | Subclass deepdrivewe, do not modify core | Self-contained example; avoids breaking changes to the reference library |
| 2026-04-14 | Collect contact maps in Phase 1 | Phase 2 CVAE needs them; avoids re-running trajectories |
| 2026-04-14 | Academy agent architecture retained as-is | All changes are in the pcoord/resampling layer, not agent topology |
| 2026-04-14 | Example placed in deepdrivewe repo | Self-contained reference for the community; demonstrates agentic development workflow |

---

## Open questions

- **αC-helix timescale:** What is the expected timescale of αC-helix motion in
  apo CDK9?  This governs how long each WE simulation segment needs to be.
  Short segments (50 ps) may not sample salt-bridge formation/breaking events.

- **Phospho-Thr186 in apo:** Confirm the biological context — active CDK9 is
  Thr186-phosphorylated in vivo.  Retain unless specifically studying the
  unphosphorylated form.

- **GBn2 loop drift:** Implicit solvent can cause unphysical loop extension in
  long runs.  Monitor with `inputs/04_verify_pcoord_residues.py` and reduce
  `rmsd_boundary_ang` if frequent recycling is observed.

- **Phase 2 contact map alignment:** Apo (321 residues) and holo (573 residues)
  produce contact maps of different sizes.  Do not train a single CVAE on both
  without alignment/padding.

- **Additional binding partners:** Are there CDK9 regulators beyond CyclinT1
  worth including in later phases (e.g., 7SK snRNP components)?
