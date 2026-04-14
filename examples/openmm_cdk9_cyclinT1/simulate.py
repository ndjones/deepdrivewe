"""CDK9 / CyclinT1 simulation components for Phase 1 weighted ensemble.

Progress coordinate layout
--------------------------
pcoord[0]  Cα RMSD to the reference PDB (Angstrom).
           Measures global structural drift; walkers that exceed
           RMSD_BOUNDARY_ANG are recycled to basis states (see recyclers.py).

pcoord[1]  Glu66 Cδ (CD) – Lys48 Nζ (NZ) Euclidean distance (Angstrom).
           Reports the αC-helix salt-bridge status:
             ~3.5 Å  → αC-in  (active / DFG-in)
             > 8 Å   → αC-out (inactive / intermediate)

Both components are computed at each reporter step and stored on
SimMetadata.pcoord as a (n_frames, 2) array per simulation segment.

Relationship to deepdrivewe
---------------------------
CDK9PcoordReporter  subclasses ContactMapRMSDReporter.
  - Inherits Cα contact-map collection (kept for Phase 2 CVAE compatibility).
  - Overrides report() to also compute the Glu66–Lys48 salt-bridge distance.
  - Overrides get_rmsds() to return shape (n_frames, 2) instead of (n_frames, 1).

CDK9SimulationAgent  subclasses SimulationAgent.
  - Overrides run_simulation() to inject CDK9PcoordReporter instead of
    ContactMapRMSDReporter.  All other logic is identical to the parent.

See CLAUDE.md for design rationale.
"""

from __future__ import annotations

import asyncio
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from academy.agent import action

from deepdrivewe import SimMetadata
from deepdrivewe.academy_agents.simulation import SimulationAgent
from deepdrivewe.simulation.openmm import ContactMapRMSDReporter
from deepdrivewe.simulation.openmm import OpenMMSimulation

try:
    import openmm
    from openmm import app
except ImportError:
    pass  # Handled at runtime


class CDK9PcoordReporter(ContactMapRMSDReporter):
    """Two-dimensional progress coordinate reporter for CDK9.

    Extends ContactMapRMSDReporter with a second pcoord dimension:
    the Glu66 CD – Lys48 NZ salt-bridge distance.  Contact maps are
    still accumulated so that CVAE analysis (Phase 2) can be added
    without re-running simulations.

    Parameters
    ----------
    report_interval : int
        Frames between reports (inherited).
    reference_file : Path
        Reference PDB for Cα RMSD superposition (inherited).
    cutoff_angstrom : float
        Contact-map cutoff in Angstrom (inherited, default 8 Å).
    mda_selection : str
        MDAnalysis atom selection for RMSD/contact map (inherited).
    openmm_selection : sequence of str
        OpenMM atom names for RMSD/contact map (inherited).
    glu_resnum : int
        Residue sequence number of Glu66 in chain A (default 66).
    lys_resnum : int
        Residue sequence number of Lys48 in chain A (default 48).
    chain_id : str
        Chain identifier for both residues (default 'A' = CDK9).
    """

    def __init__(
        self,
        report_interval: int,
        reference_file: Path,
        cutoff_angstrom: float = 8.0,
        mda_selection: str = 'protein and name CA',
        openmm_selection: tuple[str, ...] = ('CA',),
        glu_resnum: int = 66,
        lys_resnum: int = 48,
        chain_id: str = 'A',
    ) -> None:
        super().__init__(
            report_interval=report_interval,
            reference_file=reference_file,
            cutoff_angstrom=cutoff_angstrom,
            mda_selection=mda_selection,
            openmm_selection=openmm_selection,
        )
        self._glu_resnum = glu_resnum
        self._lys_resnum = lys_resnum
        self._chain_id = chain_id

        # Atom indices resolved lazily on first report() call.
        self._glu_cd_idx: int | None = None
        self._lys_nz_idx: int | None = None
        self._salt_bridge_indexed = False

        # Accumulated salt-bridge distances (one per report step).
        self._pcoord2: list[float] = []

    def _resolve_salt_bridge_indices(self, simulation: app.Simulation) -> None:
        """Identify Glu66 CD and Lys48 NZ atom indices from the topology.

        Called once on the first report() invocation.  Index lookup is
        tolerant of insertion codes in residue IDs (tries int cast).
        """
        for atom in simulation.topology.atoms():
            res = atom.residue
            try:
                resnum = int(res.id)
            except ValueError:
                continue

            chain = res.chain.id

            if chain == self._chain_id and resnum == self._glu_resnum and atom.name == 'CD':
                self._glu_cd_idx = atom.index
            elif chain == self._chain_id and resnum == self._lys_resnum and atom.name == 'NZ':
                self._lys_nz_idx = atom.index

        if self._glu_cd_idx is None:
            print(
                f'WARNING: CDK9PcoordReporter: Glu{self._glu_resnum} CD not found '
                f'in chain {self._chain_id}. pcoord[1] will be NaN.',
                flush=True,
            )
        if self._lys_nz_idx is None:
            print(
                f'WARNING: CDK9PcoordReporter: Lys{self._lys_resnum} NZ not found '
                f'in chain {self._chain_id}. pcoord[1] will be NaN.',
                flush=True,
            )
        self._salt_bridge_indexed = True

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        """Generate a two-dimensional pcoord report.

        Calls the parent to accumulate the contact map and Cα RMSD
        (pcoord[0]), then computes the salt-bridge distance (pcoord[1]).
        """
        # Parent handles contact map + RMSD accumulation.
        super().report(simulation, state)

        # Resolve atom indices once per simulation segment.
        if not self._salt_bridge_indexed:
            self._resolve_salt_bridge_indices(simulation)

        # Compute salt-bridge distance (nm → Angstrom).
        if self._glu_cd_idx is not None and self._lys_nz_idx is not None:
            positions = state.getPositions(asNumpy=True)
            glu_pos = positions[self._glu_cd_idx] * 10.0  # nm → Å
            lys_pos = positions[self._lys_nz_idx] * 10.0
            dist = float(np.linalg.norm(np.array(glu_pos) - np.array(lys_pos)))
        else:
            dist = float('nan')

        self._pcoord2.append(dist)

    def get_rmsds(self) -> np.ndarray:
        """Return the 2D progress coordinate array.

        Returns
        -------
        np.ndarray
            Shape (n_frames, 2): column 0 = Cα RMSD (Å),
            column 1 = Glu66 CD – Lys48 NZ distance (Å).
        """
        rmsd = np.array(self._rmsd)      # populated by parent report()
        sb = np.array(self._pcoord2)
        n = min(len(rmsd), len(sb))
        return np.column_stack([rmsd[:n], sb[:n]])


# ---------------------------------------------------------------------------
# CDK9 simulation agent
# ---------------------------------------------------------------------------


class CDK9SimulationAgent(SimulationAgent):
    """Simulation agent for CDK9 weighted ensemble.

    Identical to SimulationAgent except that it injects CDK9PcoordReporter
    instead of ContactMapRMSDReporter so that both pcoord dimensions are
    accumulated during each simulation segment.

    See CLAUDE.md for why we subclass rather than modify the base class.
    """

    @action
    async def run_simulation(
        self,
        metadata: dict[str, Any],
    ) -> dict[str, Any]:
        """Run an MD simulation with the CDK9 2D pcoord reporter.

        Parameters and return value are identical to
        SimulationAgent.run_simulation().  Only the reporter class is swapped.
        """
        self._log_action('run_simulation', sim_id=metadata.get('simulation_id'))

        sim_metadata = SimMetadata(**metadata)
        sim_metadata.mark_simulation_start()

        try:
            sim_output_dir = (
                self.config.output_dir / sim_metadata.simulation_name
            )

            if sim_output_dir.exists():
                await asyncio.sleep(1)
                shutil.rmtree(sim_output_dir)

            sim_output_dir.mkdir(parents=True, exist_ok=True)

            self.config.simulation_config.dump_yaml(
                sim_output_dir / 'config.yaml',
            )

            simulation = OpenMMSimulation(
                config=self.config.simulation_config,
                output_dir=sim_output_dir,
                checkpoint_file=sim_metadata.parent_restart_file,
            )

            # Inject CDK9PcoordReporter instead of ContactMapRMSDReporter.
            reporters = []
            if self.config.reference_file is not None:
                reporter = CDK9PcoordReporter(
                    report_interval=self.config.simulation_config.report_steps,
                    reference_file=self.config.reference_file,
                    cutoff_angstrom=self.config.cutoff_angstrom,
                    mda_selection=self.config.mda_selection,
                    openmm_selection=tuple(self.config.openmm_selection),
                )
                reporters.append(reporter)

            await asyncio.to_thread(simulation.run, reporters=reporters)

            if reporters:
                pcoord = reporters[0].get_rmsds()          # (n_frames, 2)
                contact_maps = reporters[0].get_contact_maps()
            else:
                pcoord = []
                contact_maps = []

            trajectory_data = {
                'restart_file': str(simulation.restart_file),
                'trajectory_file': str(simulation.trajectory_file),
                'log_file': str(simulation.log_file),
            }

            sim_metadata.restart_file = simulation.restart_file
            sim_metadata.pcoord = (
                pcoord.tolist() if hasattr(pcoord, 'tolist') else list(pcoord)
            )
            sim_metadata.mark_simulation_end()

            return {
                'metadata': sim_metadata.model_dump(),
                'trajectory': trajectory_data,
                'contact_maps': (
                    contact_maps.tolist()
                    if hasattr(contact_maps, 'tolist')
                    else list(contact_maps)
                ),
                'rmsd': (
                    pcoord.tolist()
                    if hasattr(pcoord, 'tolist')
                    else list(pcoord)
                ),
                'success': True,
            }

        except Exception as e:
            self._log_error('run_simulation', e, sim_id=metadata.get('simulation_id'))
            sim_metadata.mark_simulation_end()

            return {
                'metadata': sim_metadata.model_dump(),
                'trajectory': {},
                'contact_maps': [],
                'rmsd': [],
                'success': False,
                'error': str(e),
            }
