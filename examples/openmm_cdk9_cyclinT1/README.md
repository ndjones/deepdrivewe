# CDK9 / CyclinT1 — Weighted Ensemble Example

A deepdrivewe + Academy agents example for mapping the conformational landscape
of CDK9 and understanding how its activating partner CyclinT1 shifts that
landscape.  Built as a self-contained reference for both end users and for
AI-assisted development.

> **AI developers:** Read `CLAUDE.md` first — it covers the deepdrivewe
> architecture, CDK9 subclasses, design decisions, and how to continue
> AI-assisted work on this codebase.

---

## Scientific question

Does CyclinT1 lock CDK9's αC-helix in the active (αC-in) conformation, and if
so by how much?  Two independent weighted ensemble runs — CDK9 alone (apo) vs.
CDK9 + CyclinT1 (holo) — are compared post-hoc to answer this.

| Condition | Simulation | Expected behaviour |
|-----------|-----------|-------------------|
| `apo` | CDK9 alone (PDB 4BCI chain A) | Free αC-helix mobility; αC-in and αC-out both sampled |
| `holo_cyclinT1` | CDK9 + CyclinT1 (chains A+B) | CyclinT1 expected to stabilise αC-in |

---

## Development status

| Phase | Branch | Status | Description |
|-------|--------|--------|-------------|
| 1 | `feature/phase1-geometric-pcoord` | **Code complete** | Geometric 2D pcoord; independent apo/holo runs |
| 2 | `feature/phase2-cvae-pcoord` | Planned | CVAE latent pcoord; Voronoi binning |
| 3 | `feature/phase3-comparative-orchestration` | Planned | Coupled comparative orchestration; divergence-driven resampling |

---

## Where to look

| Need | Go to |
|------|-------|
| Run the example | [`docs/quickstart.md`](docs/quickstart.md) |
| Scientific background and goals | [`dev/intent.md`](dev/intent.md) |
| Technical design and decisions | [`dev/spec.md`](dev/spec.md) |
| Phase-specific plans and status | [`dev/phases/`](dev/phases/) |
| AI assistant context and architecture | [`CLAUDE.md`](CLAUDE.md) |
| HPC job submission | [`scripts/`](scripts/) |
| Structure preparation pipeline | [`inputs/`](inputs/) |

---

## Directory structure

```
examples/openmm_cdk9_cyclinT1/
├── README.md                      This file
├── CLAUDE.md                      AI assistant guide (read first)
├── docs/
│   └── quickstart.md              How to run: prereqs, steps, config, outputs
├── dev/
│   ├── intent.md                  Scientific goals and biological motivation
│   ├── spec.md                    Technical design decisions and decision log
│   └── phases/
│       ├── phase1.md              Phase 1 plan and deliverables (complete)
│       ├── phase2.md              Phase 2 plan: CVAE latent pcoord (planned)
│       └── phase3.md              Phase 3 plan: comparative orchestration (planned)
├── main.py                        Entry point
├── simulate.py                    CDK9PcoordReporter + CDK9SimulationAgent
├── recyclers.py                   BoundaryRecycler
├── binners.py                     Rectilinear2DBinner
├── orchestrator.py                CDK9OrchestratorAgent
├── config_apo.yaml                Apo condition config
├── config_holo_cyclinT1.yaml      Holo condition config
├── inputs/                        Structure preparation pipeline (4 scripts)
└── scripts/                       HPC Slurm job scripts
```

---

## Comparison with NTL9 Academy example

| Feature | NTL9 (`openmm_ntl9_hk_academy`) | CDK9 (`openmm_cdk9_cyclinT1`) |
|---------|--------------------------------|-------------------------------|
| pcoord dims | 1 (RMSD) | 2 (RMSD + αC-helix salt bridge) |
| Binner | 1D `RectilinearBinner` | 2D `Rectilinear2DBinner` |
| Recycler | `LowRecycler` (target-directed) | `BoundaryRecycler` (landscape mapping) |
| Reporter | `ContactMapRMSDReporter` | `CDK9PcoordReporter` (subclass) |
| Orchestrator goal | Folding progress | αC-in fraction logging |
| Scientific goal | Protein folding | Conformational landscape comparison |
