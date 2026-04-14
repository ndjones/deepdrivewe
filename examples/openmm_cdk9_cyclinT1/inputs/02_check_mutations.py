"""Inspect the three engineered mutations in CyclinT1 (chain B, PDB 4BCI).

This is step 2 of 3 in the structure preparation pipeline. Run from the inputs/
directory after 01_download_and_clean.py:

    python 02_check_mutations.py

What this script does
---------------------
CyclinT1 in 4BCI contains three mutations introduced during construct engineering
for crystallisation. This script:

1. Reports each mutation with its position and amino acid change.
2. Measures the minimum distance from each mutated residue to any CDK9 (chain A)
   residue in the holo crystal structure.
3. Classifies each mutation as interface-proximal (< 5 Å) or distal.
4. Prints a recommendation on whether to revert each mutation to wild-type.

The three mutations (from PDB SEQADV records) are:
  - Q77R:  Gln → Arg at CyclinT1 position 77
  - E96G:  Glu → Gly at CyclinT1 position 96
  - F241L: Phe → Leu at CyclinT1 position 241

Decision guidance
-----------------
If a mutation is within 5 Å of CDK9 (interface-proximal), it may alter the
CDK9/CyclinT1 interaction energy and change the conformational ensemble we are
trying to study. In that case, revert the mutation to wild-type in the holo
structure before proceeding to equilibration.

Reverting a mutation requires replacing the side chain with the canonical residue.
This can be done with:
  - PDBFixer: fixer.applyMutations({'B': ['GLN-77-GLN']}, 'B')  # Q77R → revert to Q
  - PyMOL mutagenesis wizard
  - MODELLER

If a mutation is distal (> 5 Å), it is unlikely to affect the CDK9 conformational
landscape and can be left as-is for Phase 1.

Outputs
-------
Prints a report to stdout. No files are written.

After reviewing the report, edit 03_equilibrate.py if any mutations need reverting.
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Known mutations from 4BCI SEQADV records
# ---------------------------------------------------------------------------

MUTATIONS = [
    {'position': 77,  'from_aa': 'GLN', 'to_aa': 'ARG', 'one_letter': 'Q77R'},
    {'position': 96,  'from_aa': 'GLU', 'to_aa': 'GLY', 'one_letter': 'E96G'},
    {'position': 241, 'from_aa': 'PHE', 'to_aa': 'LEU', 'one_letter': 'F241L'},
]

# Distance threshold (Å) for classifying interface proximity.
INTERFACE_CUTOFF_ANGSTROM = 5.0

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

INPUTS_DIR = Path(__file__).parent
HOLO_PDB = INPUTS_DIR / 'holo_cyclinT1' / 'cdk9_holo_cyclinT1_clean.pdb'


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------


def check_mutations(holo_pdb: Path) -> None:
    """Load the holo structure and measure mutation proximity to CDK9."""
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import distances
        import numpy as np
    except ImportError:
        print(
            'ERROR: MDAnalysis and numpy are required.\n'
            '  conda install -c conda-forge mdanalysis numpy',
            file=sys.stderr,
        )
        sys.exit(1)

    if not holo_pdb.exists():
        print(
            f'ERROR: {holo_pdb} not found.\n'
            '  Run 01_download_and_clean.py first.',
            file=sys.stderr,
        )
        sys.exit(1)

    print(f'Loading {holo_pdb.relative_to(INPUTS_DIR)} ...')
    u = mda.Universe(str(holo_pdb))

    # CDK9 is chain A; CyclinT1 is chain B.
    cdk9 = u.select_atoms('segid A')
    cyclin = u.select_atoms('segid B')

    if len(cdk9) == 0 or len(cyclin) == 0:
        # MDAnalysis may use 'chainID' instead of 'segid' depending on PDB format.
        cdk9 = u.select_atoms('chainID A')
        cyclin = u.select_atoms('chainID B')

    if len(cdk9) == 0:
        print(
            'WARNING: Could not select CDK9 (chain A). '
            'Check that the holo PDB was written with chain IDs preserved.',
        )

    print(f'  CDK9 (chain A): {len(cdk9)} atoms')
    print(f'  CyclinT1 (chain B): {len(cyclin)} atoms')

    print('\n' + '=' * 60)
    print('CyclinT1 Mutation Interface Proximity Report')
    print('=' * 60)
    print(
        f'Interface cutoff: {INTERFACE_CUTOFF_ANGSTROM} Å\n'
        f'(minimum distance from mutated residue heavy atoms to any CDK9 heavy atom)\n',
    )

    recommendations = []

    for mut in MUTATIONS:
        pos = mut['position']
        from_aa = mut['from_aa']
        to_aa = mut['to_aa']
        label = mut['one_letter']

        # Select the mutated residue in the structure (it appears as the TO_AA
        # residue since we are looking at the crystal structure as deposited).
        mut_residue = u.select_atoms(
            f'(segid B or chainID B) and resnum {pos}',
        )

        if len(mut_residue) == 0:
            print(f'  {label}: residue B{pos} not found in structure (may be missing).')
            recommendations.append((label, 'UNKNOWN', from_aa, to_aa, pos))
            continue

        # Compute minimum distance to CDK9.
        if len(cdk9) > 0:
            dist_matrix = distances.distance_array(
                mut_residue.positions,
                cdk9.positions,
            )
            min_dist = float(np.min(dist_matrix))
            at_interface = min_dist < INTERFACE_CUTOFF_ANGSTROM
        else:
            min_dist = float('nan')
            at_interface = False

        status = 'INTERFACE-PROXIMAL' if at_interface else 'distal'
        flag = ' *** REVIEW RECOMMENDED ***' if at_interface else ''

        print(f'  {label}  (CyclinT1 pos {pos}):')
        print(f'    Wild-type: {from_aa}  →  Crystal: {to_aa}')
        print(f'    Min distance to CDK9: {min_dist:.1f} Å  [{status}]{flag}')
        print()

        recommendations.append((label, status, from_aa, to_aa, pos))

    # Summary recommendation table.
    print('=' * 60)
    print('Summary')
    print('=' * 60)
    print(f'  {"Mutation":<8} {"Status":<22} {"Recommendation"}')
    print(f'  {"-"*8} {"-"*22} {"-"*40}')

    for label, status, from_aa, to_aa, pos in recommendations:
        if status == 'INTERFACE-PROXIMAL':
            rec = f'REVERT to {from_aa} — may affect CDK9 conformation'
        elif status == 'distal':
            rec = 'Accept as-is for Phase 1'
        else:
            rec = 'Manual inspection required'
        print(f'  {label:<8} {status:<22} {rec}')

    print()
    print(
        'To revert a mutation: edit the MUTATIONS_TO_REVERT list in 03_equilibrate.py\n'
        'and set it to the list of mutations you want reverted (e.g., ["Q77R"]).',
    )
    print('=' * 60)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Run the mutation proximity check."""
    print('=' * 60)
    print('CDK9 / CyclinT1 Structure Preparation — Step 2 of 3')
    print('=' * 60)
    check_mutations(HOLO_PDB)
    print('\nNext step: review report, then run 03_equilibrate.py')


if __name__ == '__main__':
    main()
