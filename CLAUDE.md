# CLAUDE.md — deepdrivewe

> **For AI coding assistants (Claude Code, Cursor, Copilot, etc.)**
> Read this file before touching any code. It provides the context needed
> to work productively on this codebase without going in circles.

## What is deepdrivewe?

DeepDriveWE is an implementation of the WESTPA weighted ensemble (WE) sampling
algorithm, extended with machine-learning-driven adaptive sampling.  It is built
on two task-distribution frameworks:

- **Colmena** (legacy): `WorkQueue`-based task server, still used in
  `examples/openmm_ntl9_hk/` and all non-`_academy` examples.
- **Academy** (active development, `feature/academy-agents` branch): async
  agent framework from Globus Labs.  All new examples should use Academy.

**This branch (`feature/academy-agents`) is the Academy branch.**
Do not use Colmena patterns here; use Academy agents.

## Architecture at a glance

```
deepdrivewe/
├── deepdrivewe/                # Python package
│   ├── api.py                  # SimMetadata, BasisStates, TargetState, pcoord types
│   ├── academy_agents/         # Academy agent classes
│   │   ├── orchestrator.py     # OrchestratorAgent — workflow loop
│   │   ├── ensemble.py         # EnsembleManagerAgent — binner/resampler/recycler
│   │   ├── simulation.py       # SimulationAgent — runs one OpenMM segment
│   │   └── config.py           # SimulationPoolConfig, AcademyWorkflowConfig
│   ├── simulation/
│   │   └── openmm.py           # OpenMMConfig, OpenMMSimulation, ContactMapRMSDReporter
│   ├── binners/
│   │   └── rectilinear.py      # RectilinearBinner (1D)
│   ├── recyclers/
│   │   └── base.py             # Recycler ABC + recycle_simulations()
│   └── resamplers/             # HuberKimResampler, etc.
└── examples/
    ├── openmm_ntl9_hk_academy/ # NTL9 folding reference example (Academy)
    └── openmm_cdk9_cyclinT1/   # CDK9/CyclinT1 Phase 1 example ← NEW
```

## How to add a new example

1. Create `examples/<name>/` — user-facing: README.md, CLAUDE.md, configs, scripts.
2. Put all Python modules in the same directory (self-contained).
3. In `main.py`, add `sys.path.insert(0, str(Path(__file__).parent))` so module
   imports (`from simulate import ...`) work when run from repo root.
4. Subclass the four extension points in `deepdrivewe/`:
   - `ContactMapRMSDReporter` → custom reporter for your pcoord
   - `SimulationAgent` → inject your reporter in `run_simulation()`
   - `Recycler` → implement `recycle(pcoords)` returning indices
   - `OrchestratorAgent` → override `evaluate_goals()` loop (optional)
5. Run from repo root: `python examples/<name>/main.py --config examples/<name>/config.yaml`

## Key subclassing patterns

### Custom progress coordinate (most common task)

```python
# In your_example/simulate.py
from deepdrivewe.simulation.openmm import ContactMapRMSDReporter

class MyReporter(ContactMapRMSDReporter):
    def report(self, simulation, state):
        super().report(simulation, state)   # keeps contact maps
        # append your extra metric to self._pcoord2 list

    def get_rmsds(self) -> np.ndarray:
        # return shape (n_frames, N) — N = number of pcoord dims
        ...
```

### Custom simulation agent

```python
from deepdrivewe.academy_agents.simulation import SimulationAgent
from academy.agent import action

class MySimulationAgent(SimulationAgent):
    @action
    async def run_simulation(self, metadata):
        # copy parent logic, swap reporter class
        ...
```

### Custom recycler

```python
from deepdrivewe.recyclers.base import Recycler

class MyRecycler(Recycler):
    def recycle(self, pcoords: np.ndarray) -> np.ndarray:
        # pcoords shape: (n_walkers, n_dims)
        # return integer indices of walkers to recycle
        return np.where(pcoords[:, 0] > self.threshold)[0]
```

## pcoord conventions

- `SimMetadata.pcoord` is a `list[list[float]]` — shape `(n_frames, n_dims)`.
- `ContactMapRMSDReporter.get_rmsds()` returns `(n_frames, 1)` by default.
- Custom reporters can return `(n_frames, N)` — the resampler/binner must be
  compatible with N dimensions.
- The last frame `pcoord[-1]` is what the recycler and binner see.

## Testing

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run full test suite
pytest tests/ -x

# Run a minimal workflow smoke test (NTL9 reference example, CPU, 3 iterations)
export OPENMM_CPU_THREADS=1
python examples/openmm_ntl9_hk_academy/main_academy.py \
    --config examples/openmm_ntl9_hk_academy/config_minimal.yaml
```

## Common mistakes to avoid

- **Don't import from `deepdrivewe.examples.*`** — the `examples/` subdirectory
  in the package is used by older Colmena examples only and is not a proper
  package (no `__init__.py`).  Use `sys.path` tricks or direct imports instead.

- **Don't copy Colmena patterns** — `PipeQueues`, `ParslTaskServer`, and
  `WESTPAThinker` are the old system.  New code uses Academy `@action`,
  `@loop`, `Handle`, and `Manager`.

- **`@action` must be on overridden methods** — if you override `run_simulation()`
  from `SimulationAgent`, keep the `@action` decorator, or the Academy framework
  won't register the method as callable via handles.

- **`get_rmsds()` shape must match binner dimension** — if your binner uses
  `pcoord_idx=1` but your reporter only returns 1 column, you'll get an
  IndexError at bin assignment time.

## Current branch status

Branch: `feature/academy-agents`
Active example: `examples/openmm_cdk9_cyclinT1/` (Phase 1 geometric pcoord)
Reference example: `examples/openmm_ntl9_hk_academy/` (NTL9 folding, Academy)

See `examples/openmm_cdk9_cyclinT1/CLAUDE.md` for CDK9-specific context
and `examples/openmm_cdk9_cyclinT1/README.md` for the full phased roadmap.
