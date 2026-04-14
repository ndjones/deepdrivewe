"""CDK9 / CyclinT1 orchestrator agent for Phase 1 weighted ensemble.

CDK9OrchestratorAgent
---------------------
Subclasses OrchestratorAgent (deepdrivewe) and overrides the
evaluate_goals loop with CDK9-specific coverage metrics.

In Phase 1 there are no goal-directed reward signals: the only objective
is uniform exploration of the (RMSD, salt-bridge) 2D space.  The
evaluate_goals loop therefore logs bin-coverage statistics at each
polling interval rather than modifying the sampling strategy.

Phase 4 extension point
-----------------------
When reward-based guidance is added, extend this class to:
1. Receive a CDK9ConformationalGoal (e.g. salt-bridge < 4 Å).
2. Query the ensemble state for fraction of probability mass in goal bins.
3. Emit a reward signal to adjust walker allocation toward goal cells.

See CLAUDE.md for the full phased roadmap.
"""

from __future__ import annotations

import asyncio

from academy.agent import loop

from deepdrivewe.academy_agents.orchestrator import OrchestratorAgent


class CDK9OrchestratorAgent(OrchestratorAgent):
    """Orchestrator with CDK9-specific goal evaluation.

    All workflow coordination (iteration loop, checkpointing, status) is
    inherited from OrchestratorAgent.  Only evaluate_goals is overridden.
    """

    @loop
    async def evaluate_goals(self, shutdown: asyncio.Event) -> None:
        """Log 2D bin-coverage statistics for the running ensemble.

        Polls the ensemble state every 60 seconds and reports:
        - Fraction of active walkers in the αC-in region (salt-bridge < 4 Å).

        This is a monitoring-only loop; it does not modify walker weights.
        Extend in Phase 4 to add adaptive guidance.

        Parameters
        ----------
        shutdown : asyncio.Event
            Event set by the Academy framework on graceful shutdown.
        """
        self.logger.info('Starting CDK9 evaluate_goals loop (coverage monitoring)')

        while not shutdown.is_set() and not self._workflow_complete:
            try:
                state = await self.ensemble_manager.get_ensemble_state()
                next_sims = state.get('next_sims', [])

                if next_sims:
                    pcoords = []
                    for sim in next_sims:
                        pcoord = sim.get('pcoord', [])
                        if pcoord:
                            last = pcoord[-1]
                            if isinstance(last, (list, tuple)) and len(last) >= 2:
                                pcoords.append((float(last[0]), float(last[1])))

                    if pcoords:
                        sb_values = [p[1] for p in pcoords]
                        ac_in_frac = sum(1 for v in sb_values if v < 4.0) / len(sb_values)
                        self.logger.info(
                            f'[CDK9] iteration={self._current_iteration} '
                            f'active_walkers={len(pcoords)} '
                            f'aC_in_fraction={ac_in_frac:.2f} (salt_bridge < 4 Å)',
                        )

                await asyncio.sleep(60.0)

            except Exception as e:
                self._log_error('evaluate_goals', e)
                await asyncio.sleep(15.0)

        self.logger.info('Exiting CDK9 evaluate_goals loop')
