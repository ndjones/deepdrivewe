# Phase 3 — Coupled Comparative Ensembles

**Branch:** `feature/phase3-comparative-orchestration`
**Status:** Planned
**Depends on:** Phase 2 CVAE latent pcoord (shared coordinate system across conditions)

## Intent

Extend the orchestrator to manage apo and holo ensembles within a single
workflow, enabling online comparison.  Adaptively focuses compute on
conformational regions where the two CDK9 landscapes diverge most — rather
than post-hoc comparison after both runs complete independently.

## 3.1 ComparativeOrchestratorAgent

New agent managing two `EnsembleManagerAgent` + `SimulationPoolAgent` pairs.
Each iteration:

1. Advances both ensembles in parallel
2. Projects both result sets into the shared CVAE latent space
3. Computes divergence between the two bin weight distributions
   (Jensen-Shannon divergence or Earth mover's distance over Voronoi bins)
4. Allocates more walkers to regions where the landscapes diverge most

## 3.2 Divergence-driven resampling

New resampler (or resampler modifier) that upweights high-divergence bins.
Walkers in those bins are split preferentially; walkers in consensus regions
are merged.  Focuses compute precisely on the conformational regions where
CyclinT1 has the largest effect on CDK9.

## 3.3 Shared CVAE

Single CVAE trained on contact maps pooled from both conditions, providing
a common latent space for comparison.

**Main technical risk:** Over-specialisation to one condition if training data
is imbalanced.  **Mitigation:** balanced mini-batches with equal sampling from
apo and holo trajectories.

**Contact map alignment:** Apo (321 residues) and holo (573 residues) produce
different-sized contact maps.  A padding or alignment strategy must be chosen
and documented before training.

## Deliverables

- [ ] `ComparativeOrchestratorAgent`
- [ ] Divergence metric computation between ensembles (JS divergence or EMD)
- [ ] Divergence-driven resampler
- [ ] Shared CVAE training pipeline with balanced apo/holo batches
- [ ] Updated analysis notebook with online divergence tracking
