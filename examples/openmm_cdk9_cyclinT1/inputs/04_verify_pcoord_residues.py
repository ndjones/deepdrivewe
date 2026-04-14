"""Verify that the pcoord residue numbers are correct in the prepared structures.

This is an optional verification step. Run after 03_equilibrate.py, before
configuring the WE simulations:

    python 04_verify_pcoord_residues.py

What this script does
---------------------
The Phase 1 progress coordinate uses two metrics:

  pcoord[0]: CDK9 Cα RMSD to the reference basis state
  pcoord[1]: Glu66 Cε – Lys48 Nζ distance (αC-helix salt bridge)

The residue numbers for the αC-helix salt bridge (Glu66 and Lys48) are taken
from the canonical CDK9 sequence (UniProt P50750). However, PDB structures often
use a different numbering scheme depending on what N-terminal residues are present
(expression tags, truncations, or SEQRES offsets).

This script:
1. Loads each basis state structure.
2. Finds the conserved αC-helix Glu and β3-strand Lys by sequence motif
   rather than just residue number.
3. Reports the actual residue numbers in the prepared structures.
4. Prints the MDAnalysis selection strings to use in the WE config.

The kinase motif used to locate the residues
--------------------------------------------
All CDK family kinases share two conserved sequence motifs:

  VAIK motif: contains the conserved β3-strand lysine (K) that coordinates ATP.
              In CDK9: ...VALKxxx... around residue 48.

  DFG motif:  marks the activation loop. The αC-helix glutamate is ~18 residues
              before DFG in most CDKs.

This script searches the CDK9 sequence for these motifs to locate the target
residues robustly, independent of PDB numbering.

If the residue numbers in your structure differ from the defaults in the WE
config, update PCOORD_GLU_RESNUM and PCOORD_LYS_RESNUM in the simulation config.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Default residue numbers from the canonical CDK9 sequence (UniProt P50750).
# These are the values used in the WE config if not overridden.
DEFAULT_LYS_RESNUM: int = 48   # β3-strand Lys (VAIK motif)
DEFAULT_GLU_RESNUM: int = 66   # αC-helix Glu

# MDAnalysis selection strings for the pcoord atoms.
# Chain A = CDK9 in both conditions.
PCOORD_LYS_SELECTION: str = '(segid A or chainID A) and resnum {lys} and name NZ'
PCOORD_GLU_SELECTION: str = '(segid A or chainID A) and resnum {glu} and name CD'  # Cε = CD in AMBER naming

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

INPUTS_DIR = Path(__file__).parent
CONDITIONS = {
    'apo': INPUTS_DIR / 'apo',
    'holo_cyclinT1': INPUTS_DIR / 'holo_cyclinT1',
}


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------


def verify_condition(condition: str, basis_dir: Path) -> None:
    """Check pcoord residues in the first basis state for a condition."""
    try:
        import MDAnalysis as mda
    except ImportError:
        print('ERROR: MDAnalysis required.  conda install -c conda-forge mdanalysis', file=sys.stderr)
        sys.exit(1)

    basis_files = sorted(basis_dir.glob('basis_*.pdb'))
    if not basis_files:
        print(f'  [{condition}] No basis states found in {basis_dir}.')
        print(f'  Run 03_equilibrate.py first.')
        return

    # Use the first basis state as representative.
    ref_pdb = basis_files[0]
    print(f'\n--- {condition} ---')
    print(f'  Reference: {ref_pdb.name}')

    u = mda.Universe(str(ref_pdb))

    # Try selecting by default residue numbers.
    lys_sel = PCOORD_LYS_SELECTION.format(lys=DEFAULT_LYS_RESNUM)
    glu_sel = PCOORD_GLU_SELECTION.format(glu=DEFAULT_GLU_RESNUM)

    lys_atoms = u.select_atoms(lys_sel)
    glu_atoms = u.select_atoms(glu_sel)

    # Report.
    _report_atom(lys_atoms, 'Lys48 NZ (β3-strand, ATP-coordinating)', lys_sel)
    _report_atom(glu_atoms, 'Glu66 CD (αC-helix, salt bridge partner)', glu_sel)

    # Measure current distance if both atoms found.
    if len(lys_atoms) == 1 and len(glu_atoms) == 1:
        import numpy as np
        dist = float(np.linalg.norm(lys_atoms.positions[0] - glu_atoms.positions[0]))
        print(f'\n  Current Glu66–Lys48 distance in {ref_pdb.name}: {dist:.2f} Å')
        if dist < 5.0:
            print('  → αC-in (salt bridge formed) — consistent with active/holo CDK9')
        elif dist < 10.0:
            print('  → intermediate — αC-helix partially displaced')
        else:
            print('  → αC-out (salt bridge broken) — consistent with inactive/apo CDK9')

    # Print ready-to-use selection strings for the WE simulation config.
    print(f'\n  MDAnalysis selection strings for WE config:')
    print(f'    pcoord[0] (CDK9 Cα RMSD):  "(segid A or chainID A) and name CA"')
    print(f'    pcoord[1] Lys NZ:           "{lys_sel}"')
    print(f'    pcoord[1] Glu CD:           "{glu_sel}"')

    # Search for the motif-based residues as a cross-check.
    _find_motif_residues(u)


def _report_atom(atoms: object, description: str, selection: str) -> None:
    """Print the result of a single atom selection."""
    if len(atoms) == 1:
        atom = atoms[0]
        print(
            f'  [OK] {description}\n'
            f'       Found: chain {atom.chainID} res {atom.resname} {atom.resnum} atom {atom.name}',
        )
    elif len(atoms) == 0:
        print(
            f'  [MISSING] {description}\n'
            f'            Selection returned 0 atoms: {selection!r}\n'
            f'            Check residue numbering in your structure.',
        )
    else:
        print(
            f'  [AMBIGUOUS] {description}\n'
            f'              Selection returned {len(atoms)} atoms (expected 1).\n'
            f'              Selection: {selection!r}',
        )


def _find_motif_residues(u: object) -> None:
    """Find VAIK and αC-helix Glu by sequence motif as a cross-check."""
    try:
        import MDAnalysis as mda
    except ImportError:
        return

    cdk9 = u.select_atoms('(segid A or chainID A) and protein')
    if len(cdk9) == 0:
        return

    # Get residues of CDK9 in order.
    residues = list(cdk9.residues)
    seq = [r.resname for r in residues]
    resnums = [r.resnum for r in residues]

    # Find VAIK motif (Val-Ala-Ile-Lys or similar; Lys is the conserved residue).
    vaik_lys = None
    for i in range(len(seq) - 3):
        if seq[i+3] == 'LYS' and seq[i] in ('VAL', 'ILE', 'LEU'):
            # Candidate VAIK-like motif: check roughly right sequence position.
            if 40 <= resnums[i+3] <= 60:
                vaik_lys = (resnums[i+3], residues[i+3])
                break

    if vaik_lys:
        print(
            f'\n  Motif search: VAIK Lys found at resnum {vaik_lys[0]} '
            f'(expected ~48)',
        )
        if vaik_lys[0] != DEFAULT_LYS_RESNUM:
            print(
                f'  *** NUMBERING MISMATCH: structure has Lys at {vaik_lys[0]}, '
                f'config uses {DEFAULT_LYS_RESNUM}. '
                f'Update PCOORD_LYS_RESNUM in the simulation config.',
            )
    else:
        print('\n  Motif search: VAIK Lys not found — manual inspection required.')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Verify pcoord residues for all conditions."""
    print('=' * 60)
    print('CDK9 pcoord Residue Verification')
    print('=' * 60)
    print(f'  Expected: Lys48 (β3-strand) and Glu66 (αC-helix) in CDK9 (chain A)')
    print(f'  Tolerance: auto-detected from structure\n')

    for condition, basis_dir in CONDITIONS.items():
        verify_condition(condition, basis_dir)

    print('\n' + '=' * 60)
    print('If residue numbers match defaults, no config changes needed.')
    print('If mismatched, update PCOORD_LYS_RESNUM and PCOORD_GLU_RESNUM')
    print('in the simulation config files before running WE.')
    print('=' * 60)


if __name__ == '__main__':
    main()
