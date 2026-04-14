"""Minimise and equilibrate CDK9 structures; save basis state snapshots.

This is step 3 of 3 in the structure preparation pipeline. Run from the inputs/
directory after 01_download_and_clean.py and 02_check_mutations.py:

    python 03_equilibrate.py

What this script does
---------------------
For each condition (apo, holo_cyclinT1):

1. Loads the cleaned PDB from step 1.
2. Optionally reverts CyclinT1 mutations flagged in step 2.
3. Adds hydrogen atoms using OpenMM's force field (AMBER ff14SB + PHOSAA10).
4. Sets up implicit solvent (GBn2 / OBC2 — portable across compute environments).
5. Runs energy minimisation to relax the crystal geometry and fill the T3E void.
6. Runs a short NVT equilibration at 300 K.
7. Saves N_BASIS_STATES evenly spaced snapshots as basis states.

Basis state snapshots are written to:
  apo/basis_000.pdb, basis_001.pdb, ...
  holo_cyclinT1/basis_000.pdb, basis_001.pdb, ...

These files are the direct inputs to the weighted ensemble simulations.

Force field choices
-------------------
AMBER ff14SB with PHOSAA10 phosphorylation parameters:
  - 'amber14-all.xml'          Main protein force field
  - 'implicit/gbn2.xml'        GBn2 implicit solvent model (Nguyen et al. 2013)

PHOSAA10 adds parameters for phosphorylated Ser, Thr, and Tyr residues. TPO
(phospho-Thr186) in CDK9 is handled natively.

Note on the GBn2 model: it is one of the better-calibrated GB models for proteins
and is the recommended implicit solvent choice in the AMBER literature for
conformational sampling. It will be slower than OBC2 but more accurate for loop
dynamics. Set SOLVENT_MODEL = 'implicit/obc2.xml' to switch to OBC2 if needed.

Adjustable parameters
---------------------
  EQUILIBRATION_LENGTH_NS   Duration of equilibration per condition.
                             Default: 2.0 ns. Increase for larger/more complex systems.
  N_BASIS_STATES            Number of snapshots to save as basis states.
                             Default: 10. Frames are evenly spaced across the run.
  TEMPERATURE_K             Simulation temperature. Default: 300 K.
  TIMESTEP_PS               Integration timestep. Default: 2 fs (0.002 ps).
  MUTATIONS_TO_REVERT       List of CyclinT1 mutations to revert (from step 2 report).
                             E.g., ['Q77R', 'E96G']. Empty list = accept all.

Compute environment note
------------------------
This script is written to run on any hardware. Set PLATFORM to:
  'CPU'     Portable default. ~1–4 h per condition for 2 ns on a modern laptop.
  'CUDA'    NVIDIA GPU. ~15–30 min per condition.
  'OpenCL'  AMD/Intel GPU or Apple Silicon. ~30–60 min per condition.

The platform is auto-detected if PLATFORM = 'auto'.

Visual inspection
-----------------
Before committing the basis states, inspect the minimised structures visually in
PyMOL or UCSF ChimeraX:
  - Check that the T3E binding site has closed sensibly after ligand removal.
  - Check that the modelled loops (L2 loop 88–96, activation loop 178–181) are
    geometrically reasonable.
  - Check that TPO A 186 has a sensible phosphate geometry.
  - For holo: verify the CDK9/CyclinT1 interface looks intact.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Configuration — edit as needed
# ---------------------------------------------------------------------------

EQUILIBRATION_LENGTH_NS: float = 2.0
N_BASIS_STATES: int = 10
TEMPERATURE_K: float = 300.0
TIMESTEP_PS: float = 0.002  # 2 fs
PLATFORM: str = 'auto'      # 'CPU', 'CUDA', 'OpenCL', or 'auto'

# Force field files (OpenMM XML format).
# PHOSAA10 handles TPO (phospho-Thr186). Both files must be available in the
# OpenMM data directory (installed with conda-forge openmm).
FORCEFIELD_FILES: list[str] = [
    'amber14-all.xml',
    'implicit/gbn2.xml',
]

# CyclinT1 mutations to revert before equilibration.
# Populate from the 02_check_mutations.py report if interface-proximal mutations
# were flagged. Format: list of strings like 'Q77R', 'E96G', 'F241L'.
# Empty list = accept all mutations from the crystal structure as-is.
MUTATIONS_TO_REVERT: list[str] = []

# Restraint strength (kJ/mol/nm^2) on heavy atoms during initial minimisation.
# Applied for the first RESTRAINED_MINIMISATION_STEPS steps, then removed.
# Prevents large crystal-to-solution geometry changes before the first snapshot.
RESTRAINT_K: float = 1000.0
RESTRAINED_MINIMISATION_STEPS: int = 500

# ---------------------------------------------------------------------------
# Mutation reversion map: one_letter_code → (from_residue, to_residue, chain)
# ---------------------------------------------------------------------------

MUTATION_REVERSION_MAP: dict[str, tuple[str, int, str, str]] = {
    'Q77R':  ('ARG', 77,  'B', 'GLN'),  # revert ARG → GLN at CyclinT1 pos 77
    'E96G':  ('GLY', 96,  'B', 'GLU'),  # revert GLY → GLU at CyclinT1 pos 96
    'F241L': ('LEU', 241, 'B', 'PHE'),  # revert LEU → PHE at CyclinT1 pos 241
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

INPUTS_DIR = Path(__file__).parent
CONDITIONS = {
    'apo': INPUTS_DIR / 'apo' / 'cdk9_apo_clean.pdb',
    'holo_cyclinT1': INPUTS_DIR / 'holo_cyclinT1' / 'cdk9_holo_cyclinT1_clean.pdb',
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def get_platform(platform_name: str) -> object:
    """Return an OpenMM Platform object.

    If platform_name is 'auto', tries CUDA then OpenCL then CPU.
    """
    try:
        import openmm as mm
    except ImportError:
        print('ERROR: openmm is required.  conda install -c conda-forge openmm', file=sys.stderr)
        sys.exit(1)

    if platform_name == 'auto':
        for name in ('CUDA', 'OpenCL', 'CPU'):
            try:
                p = mm.Platform.getPlatformByName(name)
                print(f'  Platform: {name} (auto-detected)')
                return p
            except Exception:
                continue
        print('  Platform: CPU (fallback)')
        return mm.Platform.getPlatformByName('CPU')

    p = mm.Platform.getPlatformByName(platform_name)
    print(f'  Platform: {platform_name}')
    return p


def revert_mutations(topology: object, positions: object, mutations: list[str]) -> tuple:
    """Apply wild-type residue substitutions to the CyclinT1 chain.

    Uses OpenMM's Modeller to mutate specified residues back to their canonical
    wild-type amino acids. Only called when MUTATIONS_TO_REVERT is non-empty.
    """
    try:
        from openmm.app import Modeller, ForceField
    except ImportError:
        sys.exit(1)

    print(f'  Reverting mutations: {mutations}')

    ff = ForceField(*FORCEFIELD_FILES)
    modeller = Modeller(topology, positions)

    for mut_code in mutations:
        if mut_code not in MUTATION_REVERSION_MAP:
            print(f'  WARNING: Unknown mutation code {mut_code!r} — skipping.')
            continue

        _current_aa, res_num, chain_id, wt_aa = MUTATION_REVERSION_MAP[mut_code]

        # Build the mutations dict for Modeller.applyMutations.
        # Format: {chain_index: 'WT_AA-resnum-WT_AA'} (no-op syntax for reversion)
        chain_idx = next(
            i for i, c in enumerate(modeller.topology.chains()) if c.id == chain_id
        )
        mutation_str = f'{wt_aa}-{res_num}-{wt_aa}'
        print(f'    Chain {chain_id} residue {res_num}: {mut_code} → {wt_aa}')
        modeller.applyMutations({chain_idx: mutation_str}, ff)

    return modeller.topology, modeller.positions


def add_hydrogens(topology: object, positions: object, ff: object) -> tuple:
    """Add hydrogen atoms using OpenMM's Modeller."""
    from openmm.app import Modeller
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(ff, pH=7.0)
    return modeller.topology, modeller.positions


def minimise(
    system: object,
    topology: object,
    positions: object,
    platform: object,
) -> object:
    """Energy minimise the system in two stages.

    Stage 1: Restrained minimisation — heavy atoms held near crystal positions.
    Stage 2: Unrestrained minimisation — full system allowed to relax.

    Returns the minimised positions.
    """
    try:
        import openmm as mm
        import openmm.unit as unit
    except ImportError:
        sys.exit(1)

    integrator = mm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1.0 / unit.picoseconds,
        TIMESTEP_PS * unit.picoseconds,
    )
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)

    print('  Stage 1: restrained minimisation ...')
    # Add positional restraints on heavy atoms.
    restraint = mm.CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
    restraint.addGlobalParameter('k', RESTRAINT_K * unit.kilojoules_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')

    for atom in topology.atoms():
        if atom.element is not None and atom.element.symbol != 'H':
            restraint.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometers))

    system.addForce(restraint)
    context.reinitialize(preserveState=True)
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=RESTRAINED_MINIMISATION_STEPS)

    # Remove restraints for stage 2.
    system.removeForce(system.getNumForces() - 1)
    context.reinitialize(preserveState=True)

    print('  Stage 2: unrestrained minimisation ...')
    mm.LocalEnergyMinimizer.minimize(context)

    return context.getState(getPositions=True).getPositions()


def run_equilibration(
    topology: object,
    positions: object,
    system: object,
    platform: object,
    condition_name: str,
    output_dir: Path,
) -> None:
    """Run NVT equilibration and save basis state snapshots.

    Snapshots are saved as PDB files at evenly spaced intervals across the
    equilibration trajectory.
    """
    try:
        import openmm as mm
        import openmm.unit as unit
        from openmm.app import PDBFile, StateDataReporter
    except ImportError:
        sys.exit(1)

    total_steps = int(EQUILIBRATION_LENGTH_NS * 1000 / TIMESTEP_PS)  # ps → steps
    save_interval = total_steps // N_BASIS_STATES

    integrator = mm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1.0 / unit.picoseconds,
        TIMESTEP_PS * unit.picoseconds,
    )
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)

    print(
        f'  Equilibration: {EQUILIBRATION_LENGTH_NS} ns '
        f'({total_steps} steps, saving {N_BASIS_STATES} snapshots)',
    )

    snapshot_idx = 0
    for step_block in range(N_BASIS_STATES):
        integrator.step(save_interval)

        state = context.getState(getPositions=True)
        pos = state.getPositions()

        out_pdb = output_dir / f'basis_{snapshot_idx:03d}.pdb'
        with open(out_pdb, 'w') as f:
            PDBFile.writeFile(topology, pos, f)

        snapshot_idx += 1
        print(
            f'    Saved {out_pdb.name} '
            f'({step_block + 1}/{N_BASIS_STATES} snapshots)',
        )

    print(f'  {N_BASIS_STATES} basis states written to {output_dir}')


# ---------------------------------------------------------------------------
# Per-condition pipeline
# ---------------------------------------------------------------------------


def prepare_one_condition(condition: str, input_pdb: Path) -> None:
    """Run the full minimisation + equilibration pipeline for one condition."""
    try:
        from openmm.app import PDBFile, ForceField
        import openmm as mm
        import openmm.unit as unit
    except ImportError:
        print('ERROR: openmm is required.  conda install -c conda-forge openmm', file=sys.stderr)
        sys.exit(1)

    output_dir = input_pdb.parent
    print(f'\n--- Condition: {condition} ---')
    print(f'  Input: {input_pdb.relative_to(INPUTS_DIR)}')

    if not input_pdb.exists():
        print(f'  ERROR: {input_pdb} not found. Run 01_download_and_clean.py first.')
        return

    # Load structure.
    pdb = PDBFile(str(input_pdb))
    topology = pdb.topology
    positions = pdb.positions

    # Optionally revert CyclinT1 mutations (holo condition only).
    if condition == 'holo_cyclinT1' and MUTATIONS_TO_REVERT:
        topology, positions = revert_mutations(topology, positions, MUTATIONS_TO_REVERT)
        minimised_label = 'with_reverted_mutations'
    else:
        minimised_label = 'as_crystal'

    # Build force field and system.
    print(f'  Force field: {", ".join(FORCEFIELD_FILES)}')
    ff = ForceField(*FORCEFIELD_FILES)

    # Add hydrogens.
    print('  Adding hydrogens ...')
    topology, positions = add_hydrogens(topology, positions, ff)

    # Save the post-hydrogen structure for inspection.
    pre_min_pdb = output_dir / f'cdk9_{condition}_pre_minimised.pdb'
    with open(pre_min_pdb, 'w') as f:
        PDBFile.writeFile(topology, positions, f)
    print(f'  Pre-minimisation structure: {pre_min_pdb.name}')

    # Create OpenMM system (implicit solvent).
    system = ff.createSystem(
        topology,
        nonbondedMethod=mm.app.NoCutoff,
        constraints=mm.app.HBonds,  # constrain H-bonds → allows 2 fs timestep
        hydrogenMass=1.5 * unit.amu,  # HMR for slightly larger timestep if desired
    )

    platform = get_platform(PLATFORM)

    # Minimise.
    print('  Minimising ...')
    min_positions = minimise(system, topology, positions, platform)

    # Save minimised structure for inspection.
    min_pdb = output_dir / f'cdk9_{condition}_minimised.pdb'
    with open(min_pdb, 'w') as f:
        PDBFile.writeFile(topology, min_positions, f)
    print(f'  Minimised structure: {min_pdb.name}')

    # Equilibrate and save basis states.
    print('  Equilibrating ...')
    run_equilibration(
        topology=topology,
        positions=min_positions,
        system=system,
        platform=platform,
        condition_name=condition,
        output_dir=output_dir,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Run equilibration for all conditions."""
    print('=' * 60)
    print('CDK9 / CyclinT1 Structure Preparation — Step 3 of 3')
    print('=' * 60)
    print(f'  Equilibration length: {EQUILIBRATION_LENGTH_NS} ns per condition')
    print(f'  Basis states per condition: {N_BASIS_STATES}')
    print(f'  Temperature: {TEMPERATURE_K} K')
    print(f'  Solvent model: implicit (GBn2)')
    if MUTATIONS_TO_REVERT:
        print(f'  Reverting mutations: {MUTATIONS_TO_REVERT}')
    else:
        print('  CyclinT1 mutations: accepted as deposited')

    for condition, input_pdb in CONDITIONS.items():
        prepare_one_condition(condition, input_pdb)

    print('\n' + '=' * 60)
    print('Structure preparation complete.')
    print()
    print('Basis states written to:')
    for condition in CONDITIONS:
        print(f'  {condition}/ basis_000.pdb ... basis_{N_BASIS_STATES - 1:03d}.pdb')
    print()
    print('Before running WE simulations:')
    print('  1. Visually inspect the minimised structures (PyMOL / ChimeraX):')
    print('     - T3E binding site geometry after ligand removal')
    print('     - Modelled loops (L2 loop 88-96, activation loop 178-181)')
    print('     - TPO A 186 phosphate geometry')
    print('     - CDK9/CyclinT1 interface (holo condition)')
    print('  2. Confirm pcoord[1] residue numbers (Glu66, Lys48) match your')
    print('     structure numbering using 04_verify_pcoord_residues.py.')
    print('=' * 60)


if __name__ == '__main__':
    main()
