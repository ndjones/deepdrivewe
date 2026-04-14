# Structure Preparation — CDK9 / CyclinT1 (PDB 4BCI)

This directory contains the structure preparation pipeline for Phase 1 of the
CDK9 conformational landscape project. It produces the equilibrated basis state PDB
files that are used as starting structures for weighted ensemble (WE) simulations.

---

## What is being prepared and why

**Source structure:** PDB 4BCI — CDK9 (chain A) + CyclinT1 (chain B) + inhibitor
T3E (chain C), resolved at 3.10 Å by X-ray crystallography.

Two simulation conditions are prepared:

| Condition | Chains kept | Purpose |
|---|---|---|
| **apo** | A (CDK9 only) | CDK9 without activating partner — samples full αC-helix mobility |
| **holo_cyclinT1** | A + B (CDK9 + CyclinT1) | CDK9 stabilised by CyclinT1 — expected to lock the αC-in active conformation |

The T3E inhibitor (chain C) is removed from both conditions. We are studying the
unliganded conformational landscape; the inhibitor structure provides only the protein
geometry at the active site.

The central question being addressed is whether CyclinT1 shifts the probability
distribution of CDK9 conformational states — specifically the αC-helix position
(quantified by the Glu66–Lys48 salt bridge distance in the WE simulations).

---

## Prerequisites

These scripts require the following packages, all installable via conda-forge:

```bash
conda install -c conda-forge pdbfixer openmm mdanalysis numpy
```

If you are using the project conda environment, these are already present.
Run all scripts from the `inputs/` directory.

---

## Pipeline overview

Run the three scripts in order:

```
Step 1:  01_download_and_clean.py    Download 4BCI and prepare cleaned PDB files
Step 2:  02_check_mutations.py       Inspect CyclinT1 mutations; decide whether to revert
Step 3:  03_equilibrate.py           Minimise + equilibrate; save basis state snapshots
```

After step 3 you will have:

```
inputs/
  apo/
    basis_000.pdb  ...  basis_00N.pdb   ← starting structures for CDK9 apo WE runs
  holo_cyclinT1/
    basis_000.pdb  ...  basis_00N.pdb   ← starting structures for CDK9+CyclinT1 WE runs
```

---

## Key structural decisions

### T3E ligand removal

T3E (a 2-amino-4-heteroaryl-pyrimidine CDK9 inhibitor) is removed before simulation.
We are studying the unliganded CDK9 landscape. The inhibitor-bound structure provides
the protein geometry; removing the ligand leaves a small void at the ATP binding site
that will relax during equilibration.

### Phospho-Thr186 (TPO) — retained

Thr186 in CDK9 is phosphorylated in the crystal structure (recorded as `TPO A 186`).
This is biologically correct: active CDK9 in the P-TEFb complex is Thr186-phosphorylated
by CDK7. Retaining this modification keeps the simulation consistent with the experimental
context and the activation loop geometry of the crystal structure.

The equilibration script uses AMBER ff14SB + PHOSAA10 force field parameters, which
fully support phospho-Thr. If you need to study unphosphorylated CDK9, set
`MUTATE_TPO_TO_THR = True` in `03_equilibrate.py`.

### Missing residues

CDK9 has 24 and CyclinT1 has 8 missing residues in 4BCI:

| Region | Residues | Strategy | Reason |
|---|---|---|---|
| CDK9 N-term (0–4, GPAKQ) | 5 | **Skip** | His-tag remnant; not part of CDK9 kinase domain |
| CDK9 L2 loop (88–96, KASPYNRCK) | 9 | **Model** | Loop is distant from αC-helix and pcoord atoms; adds structural completeness |
| CDK9 activation loop (178–181, KNSQ) | 4 | **Model** | Adjacent to Thr186 (TPO); important for correct activation loop geometry |
| CDK9 C-term (327–330, MLST) | 4 | **Skip** | Disordered C-terminal tail outside kinase domain |
| CyclinT1 N-term (0–7) | 8 | **Skip** | Expression construct tag remnant |

PDBFixer models missing residues by fragment insertion; modelled loops should be
inspected visually (e.g., in PyMOL or UCSF Chimera) before proceeding to equilibration.

### CyclinT1 mutations

Chain B contains 3 engineered mutations relative to the canonical human CyclinT1 sequence
(UniProt P21675). See `02_check_mutations.py` output for interface proximity analysis.
The script reports whether each mutation is within 5 Å of CDK9; if so, reverting to
wild-type is recommended before production runs.

The mutations are:
- **Q77R** — Gln→Arg at position 77
- **E96G** — Glu→Gly at position 96
- **F241L** — Phe→Leu at position 241

### Solvent model

Phase 1 uses **implicit solvent** (OpenMM GBn2 / OBC2) for portability across compute
environments. This reduces cost and removes the need for box equilibration. Implicit
solvent is adequate for mapping the conformational topology (which states exist and
their rough probabilities) but underestimates solvation effects on exposed loops.

Flag for explicit solvent in any publication-grade Phase 1 runs.

---

## Output files

Large trajectory and coordinate files are git-ignored (see `.gitignore`). Only the
final basis state PDB snapshots are committed to the repository.

```
raw/                        ← downloaded 4BCI.pdb (git-ignored)
apo/
  cdk9_apo_clean.pdb        ← PDBFixer output (git-ignored, large)
  cdk9_apo_minimised.pdb    ← post-minimisation (git-ignored)
  basis_000.pdb             ← committed basis state snapshot
  basis_001.pdb
  ...
holo_cyclinT1/
  cdk9_holo_clean.pdb       ← PDBFixer output (git-ignored)
  cdk9_holo_minimised.pdb   ← post-minimisation (git-ignored)
  basis_000.pdb             ← committed basis state snapshot
  basis_001.pdb
  ...
```
