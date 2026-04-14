"""Entry point for CDK9 / CyclinT1 Phase 1 weighted ensemble.

Usage
-----
Run from the deepdrivewe repository root:

    python examples/openmm_cdk9_cyclinT1/main.py \\
        --config examples/openmm_cdk9_cyclinT1/config_apo.yaml

    python examples/openmm_cdk9_cyclinT1/main.py \\
        --config examples/openmm_cdk9_cyclinT1/config_holo_cyclinT1.yaml

Prerequisites
-------------
1. Install deepdrivewe:  pip install -e .
2. Run the structure preparation pipeline first:
       cd examples/openmm_cdk9_cyclinT1/inputs
       python 01_download_and_clean.py
       python 03_equilibrate.py      # ~3-5 h on CPU; use HPC scripts
       python 04_verify_pcoord_residues.py
3. Edit the YAML config: set num_workers, hardware_platform, output_dir.

Agent topology
--------------
::

    CDK9OrchestratorAgent
    ├── SimulationPoolAgent
    │   ├── CDK9SimulationAgent (worker 0)
    │   ├── CDK9SimulationAgent (worker 1)
    │   └── CDK9SimulationAgent (worker N)
    └── EnsembleManagerAgent
            binner:    Rectilinear2DBinner (2D RMSD × salt-bridge grid)
            resampler: HuberKimResampler
            recycler:  BoundaryRecycler   (RMSD > threshold → recycle)

See CLAUDE.md for architecture overview and design decisions.
"""

from __future__ import annotations

import asyncio
import logging
import sys
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Allow running as a script from the repo root without installing the example
# as a package (e.g. python examples/openmm_cdk9_cyclinT1/main.py).
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import MDAnalysis
import numpy as np
from MDAnalysis.analysis import rms
from academy.exchange import LocalExchangeFactory
from academy.manager import Manager
from pydantic import Field

from deepdrivewe import BaseModel
from deepdrivewe import BasisStates
from deepdrivewe import EnsembleCheckpointer
from deepdrivewe import WeightedEnsemble
from deepdrivewe.academy_agents import AcademyWorkflowConfig
from deepdrivewe.academy_agents import EnsembleManagerAgent
from deepdrivewe.academy_agents import SimulationPoolAgent
from deepdrivewe.academy_agents.config import SimulationPoolConfig
from deepdrivewe.resamplers import HuberKimResampler
from deepdrivewe.simulation.openmm import OpenMMConfig

from binners import Rectilinear2DBinner
from orchestrator import CDK9OrchestratorAgent
from recyclers import BoundaryRecycler
from simulate import CDK9SimulationAgent


# ---------------------------------------------------------------------------
# Basis state pcoord initializer
# ---------------------------------------------------------------------------


class CDK9BasisStateInitializer(BaseModel):
    """Compute the 2D pcoord for each basis state PDB file.

    Returns [rmsd, salt_bridge_dist] so that the WeightedEnsemble assigns
    initial walkers to the correct 2D bin.

    Parameters
    ----------
    reference_file : Path
        Reference PDB for Cα RMSD superposition (typically basis_000.pdb).
    mda_selection : str
        MDAnalysis atom selection for the RMSD calculation.
    chain : str
        Chain ID containing Glu66 and Lys48 (CDK9 = 'A').
    glu66_resnum : int
        Residue sequence number of Glu66 (default 66).
    lys48_resnum : int
        Residue sequence number of Lys48 (default 48).
    """

    reference_file: Path = Field(
        description='Reference PDB for Cα RMSD superposition.',
    )
    mda_selection: str = Field(
        default='protein and name CA',
        description='MDAnalysis atom selection for RMSD.',
    )
    chain: str = Field(
        default='A',
        description='Chain ID for Glu66/Lys48 (CDK9 = A).',
    )
    glu66_resnum: int = Field(
        default=66,
        description='Residue sequence number of Glu66.',
    )
    lys48_resnum: int = Field(
        default=48,
        description='Residue sequence number of Lys48.',
    )

    def __call__(self, basis_file: str) -> list[float]:
        """Return [rmsd, salt_bridge_dist] for a single basis state.

        Parameters
        ----------
        basis_file : str
            Path to the basis state PDB file.

        Returns
        -------
        list[float]
            [pcoord[0], pcoord[1]] = [Cα RMSD (Å), Glu66 CD – Lys48 NZ (Å)].
        """
        u = MDAnalysis.Universe(basis_file)
        ref = MDAnalysis.Universe(str(self.reference_file))

        # pcoord[0]: Cα RMSD
        pos = u.select_atoms(self.mda_selection).positions
        ref_pos = ref.select_atoms(self.mda_selection).positions
        rmsd_val: float = rms.rmsd(pos, ref_pos, superposition=True)

        # pcoord[1]: Glu66 CD – Lys48 NZ salt-bridge distance
        sel_glu = (
            f'(segid {self.chain} or chainID {self.chain}) '
            f'and resnum {self.glu66_resnum} and name CD'
        )
        sel_lys = (
            f'(segid {self.chain} or chainID {self.chain}) '
            f'and resnum {self.lys48_resnum} and name NZ'
        )
        glu_atoms = u.select_atoms(sel_glu)
        lys_atoms = u.select_atoms(sel_lys)

        if len(glu_atoms) == 0 or len(lys_atoms) == 0:
            print(
                f'WARNING: CDK9BasisStateInitializer: Glu{self.glu66_resnum} CD '
                f'or Lys{self.lys48_resnum} NZ not found in {basis_file}. '
                f'Setting pcoord[1] = NaN.',
                flush=True,
            )
            salt_bridge_dist: float = float('nan')
        else:
            salt_bridge_dist = float(
                np.linalg.norm(glu_atoms.positions[0] - lys_atoms.positions[0])
            )

        return [rmsd_val, salt_bridge_dist]


# ---------------------------------------------------------------------------
# Experiment settings
# ---------------------------------------------------------------------------


class ExperimentSettings(BaseModel):
    """Full configuration for one CDK9 WE condition (apo or holo).

    Loaded from YAML via ExperimentSettings.from_yaml(path).
    """

    # --- Output ---
    output_dir: Path = Field(
        description='Root output directory for this run.',
    )

    # --- Iteration control ---
    num_iterations: int = Field(
        ge=1,
        description='Number of weighted ensemble iterations to run.',
    )
    max_retries: int = Field(
        default=2,
        ge=0,
        description='Maximum simulation retries on failure.',
    )

    # --- Basis / ensemble ---
    basis_states: BasisStates = Field(
        description='Basis state configuration (dir, extension, count).',
    )
    basis_state_initializer: CDK9BasisStateInitializer = Field(
        description='Computes initial pcoord for each basis state.',
    )

    # --- Binner ---
    rmsd_bin_edges: list[float] = Field(
        default=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, float('inf')],
        description='Bin edges for pcoord[0] (Cα RMSD, Å).',
    )
    salt_bridge_bin_edges: list[float] = Field(
        default=[
            0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
            float('inf'),
        ],
        description='Bin edges for pcoord[1] (salt-bridge distance, Å).',
    )
    sims_per_bin: int = Field(
        default=4,
        ge=1,
        description='Target walker count per 2D bin cell.',
    )

    # --- Resampler ---
    max_allowed_weight: float = Field(
        default=0.25,
        description='Maximum walker weight before splitting.',
    )
    min_allowed_weight: float = Field(
        default=1e-40,
        description='Minimum walker weight before merging.',
    )

    # --- Recycler ---
    rmsd_boundary_ang: float = Field(
        default=6.0,
        description='Cα RMSD threshold (Å) above which walkers are recycled.',
    )

    # --- OpenMM simulation ---
    simulation_config: OpenMMConfig = Field(
        description='OpenMM simulation parameters.',
    )
    reference_file: Path = Field(
        description='Reference PDB for RMSD / contact map calculation.',
    )
    cutoff_angstrom: float = Field(
        default=8.0,
        description='Contact map cutoff distance (Å).',
    )
    mda_selection: str = Field(
        default='protein and name CA',
        description='MDAnalysis selection for RMSD / contact map.',
    )
    openmm_selection: list[str] = Field(
        default=['CA'],
        description='OpenMM atom names for RMSD / contact map.',
    )

    # --- Academy workers ---
    num_workers: int = Field(
        default=4,
        ge=1,
        description='Number of parallel CDK9SimulationAgent workers.',
    )


# ---------------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------------


async def run_academy_workflow(cfg: ExperimentSettings) -> None:
    """Launch Academy agents and run the CDK9 weighted ensemble."""
    logging.info('Starting CDK9 Phase 1 weighted ensemble workflow')

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    checkpointer = EnsembleCheckpointer(output_dir=cfg.output_dir)
    checkpoint = checkpointer.latest_checkpoint()

    if checkpoint is None:
        ensemble = WeightedEnsemble(
            basis_states=cfg.basis_states,
            target_states=[],   # No target state: BoundaryRecycler handles boundaries
        )
        ensemble.initialize_basis_states(cfg.basis_state_initializer)
        logging.info('Initialised new weighted ensemble from basis states')
    else:
        ensemble = checkpointer.load(checkpoint)
        logging.info(f'Resumed ensemble from checkpoint: {checkpoint}')

    logging.info(f'Ensemble size: {len(ensemble.next_sims)} walkers')

    binner = Rectilinear2DBinner(
        bins_dim0=cfg.rmsd_bin_edges,
        bins_dim1=cfg.salt_bridge_bin_edges,
        bin_target_counts=cfg.sims_per_bin,
    )

    resampler = HuberKimResampler(
        sims_per_bin=cfg.sims_per_bin,
        max_allowed_weight=cfg.max_allowed_weight,
        min_allowed_weight=cfg.min_allowed_weight,
    )

    recycler = BoundaryRecycler(
        basis_states=ensemble.basis_states,
        rmsd_threshold=cfg.rmsd_boundary_ang,
    )

    sim_pool_config = SimulationPoolConfig(
        num_workers=cfg.num_workers,
        max_retries=cfg.max_retries,
        retry_delay=1.0,
        output_dir=cfg.output_dir / 'simulations',
        simulation_config=cfg.simulation_config,
        reference_file=cfg.reference_file,
        cutoff_angstrom=cfg.cutoff_angstrom,
        mda_selection=cfg.mda_selection,
        openmm_selection=cfg.openmm_selection,
    )

    workflow_config = AcademyWorkflowConfig(
        output_dir=cfg.output_dir,
        num_iterations=cfg.num_iterations,
        checkpoint_interval=1,
        simulation_pool_config=sim_pool_config,
    )

    # num_workers + 3 agents (pool + ensemble manager + orchestrator)
    num_executor_workers = cfg.num_workers + 3

    async with await Manager.from_exchange_factory(
        factory=LocalExchangeFactory(),
        executors=ThreadPoolExecutor(max_workers=num_executor_workers),
    ) as manager:
        logging.info('Academy Manager started')

        workers = []
        for i in range(cfg.num_workers):
            worker = await manager.launch(
                CDK9SimulationAgent,
                kwargs={'config': sim_pool_config},
            )
            workers.append(worker)
            logging.info(f'Launched CDK9SimulationAgent worker {i}')

        simulation_pool = await manager.launch(
            SimulationPoolAgent,
            kwargs={'config': sim_pool_config, 'workers': workers},
        )
        logging.info('Launched SimulationPoolAgent')

        ensemble_manager = await manager.launch(
            EnsembleManagerAgent,
            kwargs={
                'ensemble': ensemble,
                'binner': binner,
                'resampler': resampler,
                'recycler': recycler,
            },
        )
        logging.info('Launched EnsembleManagerAgent')

        orchestrator = await manager.launch(
            CDK9OrchestratorAgent,
            kwargs={
                'config': workflow_config,
                'simulation_pool': simulation_pool,
                'ensemble_manager': ensemble_manager,
                'checkpointer': checkpointer,
            },
        )
        logging.info('Launched CDK9OrchestratorAgent')

        await orchestrator.start_workflow()

        for iteration in range(cfg.num_iterations):
            logging.info(f'Iteration {iteration + 1}/{cfg.num_iterations}')
            success = await orchestrator.advance_iteration()
            if not success:
                logging.info('Workflow signalled complete')
                break
            status = await orchestrator.get_status()
            logging.info(f'Status: {status}')

        logging.info('All agents shutting down')

    logging.info('CDK9 Phase 1 workflow complete.')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Parse config and run the workflow."""
    parser = ArgumentParser(
        description='CDK9/CyclinT1 Phase 1 weighted ensemble (geometric pcoord)',
    )
    parser.add_argument(
        '-c',
        '--config',
        required=True,
        help='Path to YAML config (config_apo.yaml or config_holo_cyclinT1.yaml)',
    )
    args = parser.parse_args()

    cfg = ExperimentSettings.from_yaml(args.config)

    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    cfg.dump_yaml(cfg.output_dir / 'params.yaml')

    logging.basicConfig(
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        level=logging.INFO,
        handlers=[
            logging.FileHandler(cfg.output_dir / 'runtime.log'),
            logging.StreamHandler(sys.stdout),
        ],
    )

    logging.info('=' * 72)
    logging.info('CDK9 / CyclinT1 Phase 1 — Geometric pcoord WE')
    logging.info('=' * 72)
    logging.info(f'Config:        {args.config}')
    logging.info(f'Output:        {cfg.output_dir}')
    logging.info(f'Iterations:    {cfg.num_iterations}')
    logging.info(f'Workers:       {cfg.num_workers}')
    logging.info(f'RMSD boundary: {cfg.rmsd_boundary_ang} Å')
    logging.info('=' * 72)

    try:
        asyncio.run(run_academy_workflow(cfg))
        logging.info('Workflow completed successfully.')
    except Exception as e:
        logging.error(f'Workflow failed: {e}', exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
