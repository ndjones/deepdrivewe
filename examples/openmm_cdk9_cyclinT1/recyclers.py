"""Recycling logic for CDK9 / CyclinT1 weighted ensemble (Phase 1).

BoundaryRecycler
----------------
Recycles walkers whose Cα RMSD (pcoord[0]) exceeds a configurable
threshold.  This prevents trajectories from drifting into unphysically
unfolded states while keeping the WE sampling unbiased with respect to
the αC-helix salt-bridge coordinate (pcoord[1]).

There is NO target-directed recycling in Phase 1: walkers are never
preferentially restarted from a specific conformation.  Recycled walkers
are restarted from a randomly selected basis state (handled by the base
class Recycler.recycle_simulations()).

Design rationale
----------------
A threshold of 6 Å was chosen as a conservative upper bound.  The
CDK9 apo structure spans roughly 0–5 Å RMSD during normal αC-helix
transitions; values above 6 Å indicate global unfolding rather than
biologically relevant conformational change.

See CLAUDE.md for why BoundaryRecycler is used instead of LowRecycler.
"""

from __future__ import annotations

import numpy as np

from deepdrivewe.api import BasisStates
from deepdrivewe.recyclers.base import Recycler


class BoundaryRecycler(Recycler):
    """Recycle walkers that exceed the RMSD boundary.

    Parameters
    ----------
    basis_states : BasisStates
        Basis states used to restart recycled walkers.
    rmsd_threshold : float
        Cα RMSD threshold in Angstrom (default 6.0 Å).
        Walkers with pcoord[0] > threshold are recycled.
    """

    def __init__(
        self,
        basis_states: BasisStates,
        rmsd_threshold: float = 6.0,
    ) -> None:
        super().__init__(basis_states=basis_states)
        self.rmsd_threshold = rmsd_threshold

    def recycle(self, pcoords: np.ndarray) -> np.ndarray:
        """Return indices of walkers to recycle.

        Parameters
        ----------
        pcoords : np.ndarray
            Progress coordinates of shape (n_walkers, n_dims).
            Column 0 must be the Cα RMSD in Angstrom.

        Returns
        -------
        np.ndarray
            Integer indices of walkers whose RMSD exceeds the threshold.
        """
        return np.where(pcoords[:, 0] > self.rmsd_threshold)[0]
