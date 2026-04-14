# Scientific Intent

## Goal

Map the conformational landscape of CDK9 and understand how its obligate
activating partner CyclinT1 shifts that landscape, with the aim of explaining
drug selectivity mechanisms.

## Biological motivation

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

## System: PDB 4BCI

**Structure:** CDK9 (chain A, 331 residues) + CyclinT1 (chain B, 260 residues)
+ inhibitor T3E (chain C) from *Homo sapiens*, resolved at 3.10 Å by X-ray
diffraction.

**For simulation:**
- **Remove:** T3E ligand (studying the unliganded landscape)
- **Holo condition:** CDK9 (chain A) + CyclinT1 (chain B)
- **Apo condition:** CDK9 (chain A) alone

| Feature | Detail | Action |
|---------|--------|--------|
| TPO (phosphoThr186) | Activation loop, phosphorylated in crystal | Retain — reflects active CDK9; requires AMBER ff14SB + PHOSAA10 |
| CyclinT1 mutations | Q77R/E96G/F241L in chain B (crystallisation construct) | All distal (13–29 Å from CDK9 interface); accepted as-is for Phase 1 |
| Resolution 3.10 Å | Moderate; some loops poorly resolved | PDBFixer adds missing atoms/hydrogens |
| Duplicate chain IDs | PDB has A/B/A/B — second pair is HETATM | Index-based removal in `inputs/01_download_and_clean.py` |

## What success looks like

- Statistically sampled probability weight distributions over CDK9 conformational
  space for both apo and holo conditions
- Quantitative comparison showing whether CyclinT1 shifts the αC-helix
  equilibrium toward the active (αC-in) state, and by how much
- Contact map data collected and stored for Phase 2 CVAE training
- A reproducible, portable pipeline that can be re-run or extended to other
  CDK family members or binding partners

## Open questions

- **αC-helix timescale:** What is the expected timescale of αC-helix motion in
  apo CDK9?  Short WE segments (50 ps) may not capture salt-bridge
  formation/breaking events; longer segments or more iterations may be required.

- **Phospho-Thr186:** Retain phosphorylation for active CDK9 — revisit only if
  studying the unphosphorylated form specifically.

- **GBn2 loop drift:** Implicit solvent can cause unphysical loop extension in
  long runs.  Monitor recycling frequency; reduce `rmsd_boundary_ang` to 4 Å
  if walkers are recycled too often in early iterations.

- **Phase 2 contact map alignment:** Apo (321 residues) and holo (573 residues)
  produce different-sized contact maps.  Alignment or padding is required before
  training a joint CVAE on both conditions.

- **Additional binding partners:** Are CDK9 regulators beyond CyclinT1 (e.g.,
  7SK snRNP components) worth including in later phases?
