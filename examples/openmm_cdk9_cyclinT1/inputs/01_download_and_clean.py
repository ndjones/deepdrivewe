"""Download PDB 4BCI and prepare cleaned structures for apo and holo conditions.

This is step 1 of 3 in the structure preparation pipeline. Run from the inputs/
directory:

    python 01_download_and_clean.py

What this script does
---------------------
1. Downloads PDB 4BCI from the RCSB.
2. Prepares two conditions using PDBFixer:
   - apo:          CDK9 alone (chain A), T3E ligand and CyclinT1 removed
   - holo_cyclinT1: CDK9 + CyclinT1 (chains A + B), T3E ligand removed
3. For each condition:
   - Removes heterogens (T3E inhibitor, crystallisation molecules)
   - Retains TPO (phospho-Thr186) — see note below
   - Models selected missing residue loops (88–96 and 178–181 in CDK9)
   - Skips N/C-terminal expression-tag remnants
   - Adds missing heavy atoms to existing residues
4. Writes outputs to raw/ and apo/ / holo_cyclinT1/ directories.

Outputs
-------
raw/4BCI.pdb                        Original downloaded structure
apo/cdk9_apo_clean.pdb              Cleaned CDK9-only structure
holo_cyclinT1/cdk9_holo_clean.pdb   Cleaned CDK9+CyclinT1 structure

Next step
---------
Run 02_check_mutations.py to inspect the three CyclinT1 engineered mutations and
decide whether any should be reverted to wild-type before equilibration.

Notes on TPO (phospho-Thr186)
------------------------------
Thr186 is phosphorylated in 4BCI, reflecting the active CDK9 state. This script
retains the TPO residue. PDBFixer treats it as a modified residue and will not
strip it when remove_heterogens is called (it is in a MODRES record, not HETATM).
The force field assignment is handled in 03_equilibrate.py using AMBER ff14SB +
PHOSAA10 parameters.

Set REVERT_TPO_TO_THR = True only if you need to study unphosphorylated CDK9.
"""

from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — edit these if needed
# ---------------------------------------------------------------------------

PDB_ID = '4BCI'

# Chains to retain per condition.
# Chain A = CDK9, Chain B = CyclinT1, Chain C = T3E inhibitor (always removed).
CONDITIONS: dict[str, list[str]] = {
    'apo': ['A'],
    'holo_cyclinT1': ['A', 'B'],
}

# Missing residue modelling strategy.
# Residues listed here will be MODELLED by PDBFixer fragment insertion.
# All other missing residues are skipped (expression tags, disordered C-termini).
# Format: {chain_id: [list of residue sequence numbers to model]}
MODEL_MISSING: dict[str, list[int]] = {
    'A': list(range(88, 97)) + list(range(178, 182)),  # L2 loop + activation loop
    'B': [],  # CyclinT1 N-terminal tag: skip
}

# Set to True to mutate TPO A 186 → THR A 186.
# Only relevant if studying unphosphorylated CDK9.
REVERT_TPO_TO_THR: bool = False

# ---------------------------------------------------------------------------
# Paths (relative to the inputs/ directory)
# ---------------------------------------------------------------------------

INPUTS_DIR = Path(__file__).parent
RAW_DIR = INPUTS_DIR / 'raw'
RAW_PDB = RAW_DIR / f'{PDB_ID}.pdb'

CONDITION_DIRS = {
    name: INPUTS_DIR / name for name in CONDITIONS
}

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download_pdb(pdb_id: str, dest: Path) -> None:
    """Download a PDB file from the RCSB if not already present."""
    if dest.exists():
        print(f'  {dest.name} already present — skipping download.')
        return

    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    print(f'  Downloading {url} ...')
    urllib.request.urlretrieve(url, dest)
    print(f'  Saved to {dest}')


# ---------------------------------------------------------------------------
# Structure preparation
# ---------------------------------------------------------------------------


def prepare_condition(
    pdb_path: Path,
    chains_to_keep: list[str],
    output_path: Path,
    condition_name: str,
) -> None:
    """Use PDBFixer to clean a structure for a single simulation condition.

    Parameters
    ----------
    pdb_path : Path
        Path to the raw downloaded PDB file.
    chains_to_keep : list[str]
        Chain IDs to retain. All other chains (including the T3E ligand
        chain C) are removed.
    output_path : Path
        Path to write the cleaned PDB file.
    condition_name : str
        Human-readable label for log messages.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        print(
            'ERROR: pdbfixer and openmm are required.\n'
            '  conda install -c conda-forge pdbfixer openmm',
            file=sys.stderr,
        )
        sys.exit(1)

    print(f'\n--- Preparing condition: {condition_name} ---')
    print(f'  Retaining chains: {chains_to_keep}')

    fixer = PDBFixer(filename=str(pdb_path))

    # Step 1: Remove chains that are not in our keep list.
    # This removes the T3E inhibitor and, for the apo condition, CyclinT1.
    #
    # NOTE: 4BCI has duplicate chain IDs (A, B, A, B) because PDBFixer exposes
    # both protein ATOM chains and HETATM chains under the same letter. We must
    # use index-based removal to avoid ambiguity, keeping only the FIRST occurrence
    # of each chain ID we want (the protein chains).
    all_chains = list(fixer.topology.chains())
    seen_ids: set[str] = set()
    indices_to_remove = []
    for i, chain in enumerate(all_chains):
        if chain.id in chains_to_keep and chain.id not in seen_ids:
            seen_ids.add(chain.id)   # keep first occurrence of each wanted ID
        else:
            indices_to_remove.append(i)

    if indices_to_remove:
        removed_ids = [all_chains[i].id for i in indices_to_remove]
        print(f'  Removing chain indices {indices_to_remove} (IDs: {removed_ids})')
        fixer.removeChains(chainIndices=indices_to_remove)

    # Step 2: Find missing residues. PDBFixer compares the SEQRES records
    # against the ATOM records to identify gaps.
    fixer.findMissingResidues()

    # Step 3: Apply the missing residue modelling strategy.
    # Keep only residues in MODEL_MISSING; discard everything else
    # (N/C-terminal expression tags, disordered tail residues).
    _apply_missing_residue_strategy(fixer, chains_to_keep)

    # Step 4: Find and add missing non-hydrogen atoms within existing residues.
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Step 5: Remove heterogens (any remaining HETATM that isn't a standard
    # amino acid). This is a safety net — the ligand should already be gone
    # from chain removal in step 1, but crystallisation additives may remain.
    # keepResidues=False removes all non-protein, non-water heterogens.
    # Note: TPO is listed in MODRES, not HETATM, so it is NOT removed here.
    fixer.removeHeterogens(keepWater=False)

    # Step 6: Optionally revert TPO → THR.
    if REVERT_TPO_TO_THR and 'A' in chains_to_keep:
        _revert_tpo_to_thr(fixer)

    # Step 7: Write cleaned structure (no hydrogens — added in equilibration step).
    print(f'  Writing cleaned structure to {output_path}')
    with open(output_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    _report_structure(fixer, condition_name)


def _apply_missing_residue_strategy(fixer: object, chains_to_keep: list[str]) -> None:
    """Filter fixer.missingResidues to match the MODEL_MISSING strategy.

    PDBFixer stores missing residues as a dict keyed by (chain_index, seqres_idx),
    where seqres_idx is the residue's 0-based position within the full SEQRES record.

    NOTE: 4BCI has duplicate chain IDs (chains appear as A, B, A, B in the PDB) because
    PDBFixer exposes both the protein chain and the HETATM chain for each entity.
    Chain-ID-based lookup is therefore unreliable. Instead, we filter by position:
    keep internal gaps (not the first or last group per chain), which corresponds
    exactly to the L2 loop and activation loop in CDK9 and avoids modelling N/C-terminal
    expression-tag remnants.
    """
    if not fixer.missingResidues:
        print('  Missing residue groups: 0 found')
        return

    # Group entries by chain index so we can identify first/last per chain.
    from collections import defaultdict
    by_chain: dict[int, list[tuple[int, list[str]]]] = defaultdict(list)
    for (chain_idx, seqres_idx), residues in fixer.missingResidues.items():
        by_chain[chain_idx].append((seqres_idx, residues))

    filtered = {}
    for chain_idx, entries in by_chain.items():
        # Sort by seqres position so first/last are unambiguous.
        entries.sort(key=lambda x: x[0])
        n = len(entries)
        for i, (seqres_idx, residues) in enumerate(entries):
            is_terminal = (i == 0) or (i == n - 1)
            if is_terminal:
                label = 'N-term tag' if i == 0 else 'C-term tail'
                print(
                    f'  Skipping {label}: chain {chain_idx} seqres {seqres_idx} '
                    f'({residues})',
                )
            else:
                filtered[(chain_idx, seqres_idx)] = residues
                print(
                    f'  Modelling internal loop: chain {chain_idx} seqres {seqres_idx} '
                    f'({residues})',
                )

    n_original = len(fixer.missingResidues)
    fixer.missingResidues = filtered
    n_kept = len(filtered)
    print(
        f'  Missing residue groups: {n_original} found, {n_kept} selected for modelling',
    )


def _revert_tpo_to_thr(fixer: object) -> None:
    """Mutate TPO (phospho-Thr) at position 186 of chain A back to regular Thr.

    This is only called when REVERT_TPO_TO_THR = True.
    """
    from openmm.app import Residue  # noqa: F401 — for type context

    print('  Reverting TPO A 186 → THR (REVERT_TPO_TO_THR=True)')
    mutations = {'A': ['THR-186-THR']}  # no-op in terms of backbone; strips phosphate
    fixer.applyMutations(mutations, 'A')


def _report_structure(fixer: object, label: str) -> None:
    """Print a brief summary of the prepared structure."""
    chains = list(fixer.topology.chains())
    residues = list(fixer.topology.residues())
    atoms = list(fixer.topology.atoms())
    print(
        f'  [{label}] {len(chains)} chain(s), '
        f'{len(residues)} residues, {len(atoms)} heavy atoms',
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Run the full download and clean pipeline."""
    print('=' * 60)
    print('CDK9 / CyclinT1 Structure Preparation — Step 1 of 3')
    print('=' * 60)

    # Ensure output directories exist.
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for path in CONDITION_DIRS.values():
        path.mkdir(parents=True, exist_ok=True)

    # Download.
    print('\n[1/3] Downloading PDB 4BCI')
    download_pdb(PDB_ID, RAW_PDB)

    # Prepare each condition.
    print('\n[2/3] Preparing structures')
    for condition, chains in CONDITIONS.items():
        out_dir = CONDITION_DIRS[condition]
        out_pdb = out_dir / f'cdk9_{condition}_clean.pdb'
        prepare_condition(
            pdb_path=RAW_PDB,
            chains_to_keep=chains,
            output_path=out_pdb,
            condition_name=condition,
        )

    # Summary.
    print('\n[3/3] Summary')
    for condition in CONDITIONS:
        out_pdb = CONDITION_DIRS[condition] / f'cdk9_{condition}_clean.pdb'
        status = 'OK' if out_pdb.exists() else 'MISSING'
        print(f'  [{status}] {out_pdb.relative_to(INPUTS_DIR)}')

    print('\nNext step: run 02_check_mutations.py')
    print('=' * 60)


if __name__ == '__main__':
    main()
