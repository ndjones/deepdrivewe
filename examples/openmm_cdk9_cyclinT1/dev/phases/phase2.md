# Phase 2 — CVAE Latent Pcoord, Adaptive Binning

**Branch:** `feature/phase2-cvae-pcoord`
**Status:** Planned
**Depends on:** Phase 1 trajectory data (contact maps) for CVAE warm-start

## Intent

Replace the hand-crafted 2D geometric pcoord with CVAE latent coordinates
learned from contact maps collected in Phase 1.  Removes the need to
pre-specify which geometric features matter; the pcoord adapts to the actual
variance in CDK9 conformational data.

The `TrainingAgent` and `InferenceAgent` (present in deepdrivewe but unused in
Phase 1) are built for exactly this role.

## Contact map definition for CDK9

- Intra-CDK9 Cα contacts only — keeps feature space consistent between apo
  and holo conditions (CDK9 Cα count is the same in both)
- Optionally add CDK9–CyclinT1 interface contacts as a second channel in the
  holo condition
- `cutoff_angstrom`: 8.0 Å (existing deepdrivewe default)
- Contact maps are already being collected during Phase 1 runs and stored in
  `SimResult.data['contact_maps']`

**Note:** Apo (321 residues) and holo (573 residues) produce different-sized
contact maps.  Alignment or padding is required before training a joint CVAE
on both conditions.

## 2.1 CVAE warm-start

Pre-train the CVAE on Phase 1 contact maps before beginning Phase 2 WE runs.
Avoids the cold-start problem where early iterations produce poor latent
embeddings and the binner cannot effectively distribute walkers.

- Training script: `analysis/train_cvae_warmstart.py` (to be created)
- Input: contact maps from `runs/cdk9-apo/` and `runs/cdk9-holo/`
- Output: saved CVAE weights for Phase 2 initialisation

## 2.2 Enable TrainingAgent + InferenceAgent

Wire `training_agent_config` and `inference_agent_config` into
`AcademyWorkflowConfig` for the CDK9 example.  The `TrainingAgent` receives
contact maps from `SimulationAgent` and updates the CVAE online; the
`InferenceAgent` projects new snapshots into latent space and provides
latent coordinates as pcoord.

## 2.3 Voronoi binning

Replace `Rectilinear2DBinner` with a `VoronoiBinner` (new class to be added to
`deepdrivewe/binners/`) that places bins as Voronoi cells centred on previously
visited latent points.  Naturally adapts to the explored region without
requiring pre-specified grid edges.

## Deliverables

- [ ] CVAE warm-start training script (`analysis/train_cvae_warmstart.py`)
- [ ] `TrainingAgent` + `InferenceAgent` wired into CDK9 workflow
- [ ] `VoronoiBinner` in `deepdrivewe/binners/`
- [ ] Updated `config_apo.yaml` + `config_holo_cyclinT1.yaml` for CVAE pcoord
- [ ] Analysis notebook: latent space comparison apo vs. holo
