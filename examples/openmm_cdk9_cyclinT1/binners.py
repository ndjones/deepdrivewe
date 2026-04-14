"""2D rectilinear binner for CDK9 / CyclinT1 weighted ensemble (Phase 1).

Rectilinear2DBinner
-------------------
Creates a uniform 2D grid over (pcoord[0], pcoord[1]) and maps each
walker to a flat bin index:

    flat_index = row * n_cols + col

where row = np.digitize(pcoord[0], bins_dim0) and
      col = np.digitize(pcoord[1], bins_dim1).

This extends RectilinearBinner (deepdrivewe) from 1D to 2D without
modifying the upstream library.  In Phase 2 this class will be replaced
by a Voronoi binner over the CVAE latent space.

Default grid edges
------------------
pcoord[0] (Cα RMSD, Å):          [0, 1, 2, 3, 4, 5, 6, inf]   → 7 bins
pcoord[1] (salt-bridge dist, Å): [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, inf] → 11 bins
Total: 77 cells

See CLAUDE.md for design rationale.
"""

from __future__ import annotations

import numpy as np

from deepdrivewe.binners.base import Binner


class Rectilinear2DBinner(Binner):
    """2D rectilinear binner over pcoord[0] and pcoord[1].

    Parameters
    ----------
    bins_dim0 : list[float]
        Bin edges for pcoord[0] (Cα RMSD, Å).
        Must be strictly increasing; include ``float('inf')`` as last edge.
    bins_dim1 : list[float]
        Bin edges for pcoord[1] (salt-bridge distance, Å).
        Same requirements as bins_dim0.
    bin_target_counts : int
        Target number of walkers per 2D cell.
    """

    def __init__(
        self,
        bins_dim0: list[float],
        bins_dim1: list[float],
        bin_target_counts: int,
    ) -> None:
        if not np.all(np.diff(bins_dim0) > 0):
            raise ValueError('bins_dim0 must be strictly increasing.')
        if not np.all(np.diff(bins_dim1) > 0):
            raise ValueError('bins_dim1 must be strictly increasing.')

        self.bins_dim0 = bins_dim0
        self.bins_dim1 = bins_dim1
        self._bin_target_counts = bin_target_counts

    @property
    def n_bins_dim0(self) -> int:
        """Number of bins along pcoord[0]."""
        return len(self.bins_dim0) - 1

    @property
    def n_bins_dim1(self) -> int:
        """Number of bins along pcoord[1]."""
        return len(self.bins_dim1) - 1

    @property
    def nbins(self) -> int:
        """Total number of 2D cells."""
        return self.n_bins_dim0 * self.n_bins_dim1

    def get_bin_target_counts(self) -> list[int]:
        """Return per-cell target walker counts (uniform)."""
        return [self._bin_target_counts] * self.nbins

    def assign_bins(self, pcoords: np.ndarray) -> np.ndarray:
        """Map each walker to a flat 2D bin index.

        Parameters
        ----------
        pcoords : np.ndarray
            Shape (n_walkers, n_dims).  Columns 0 and 1 are used.

        Returns
        -------
        np.ndarray
            Integer flat bin indices, shape (n_walkers,).
        """
        row = np.digitize(pcoords[:, 0], self.bins_dim0, right=True)
        col = np.digitize(pcoords[:, 1], self.bins_dim1, right=True)

        # Clamp to valid range to guard against out-of-bounds values.
        row = np.clip(row, 0, self.n_bins_dim0 - 1)
        col = np.clip(col, 0, self.n_bins_dim1 - 1)

        return row * self.n_bins_dim1 + col
