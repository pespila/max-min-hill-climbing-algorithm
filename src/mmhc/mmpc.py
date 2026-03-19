"""MMPC (Max-Min Parents and Children) skeleton learning.

Replaces the C++ Forward, Backward, MaxMinHeuristic, and mmpc functions.
Uses itertools.combinations for conditioning subsets instead of UpdateCPC.
"""

from __future__ import annotations

from itertools import combinations
from typing import TYPE_CHECKING

from .independence import g_test

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


class MMPC:
    """Max-Min Parents and Children algorithm for skeleton discovery.

    Identifies the undirected skeleton of a Bayesian network by testing
    conditional independence between all variable pairs.
    """

    def __init__(
        self,
        data: NDArray[np.int_],
        cardinalities: NDArray[np.int_],
        alpha: float = 0.05,
    ) -> None:
        self._data = data
        self._cardinalities = cardinalities
        self._alpha = alpha
        self._n_vars = data.shape[1]

    def run(self) -> list[list[int]]:
        """Run MMPC and return symmetrized parent-children sets (0-indexed)."""
        pc_sets: list[list[int]] = []
        for t in range(self._n_vars):
            cpc = self._forward(t)
            cpc = self._backward(t, cpc)
            pc_sets.append(cpc)

        # Symmetry constraint: x in PC(t) only if t in PC(x)
        pc_sets = self._symmetrize(pc_sets)
        return pc_sets

    def _forward(self, target: int) -> list[int]:
        """Forward phase: iteratively add the most dependent variable to CPC."""
        cpc: list[int] = []
        remaining = [v for v in range(self._n_vars) if v != target]

        while remaining:
            best_var = -1
            best_min_assoc = -1.0
            best_g_stat = -1.0

            for x in remaining:
                min_pval, max_g = self._max_min_heuristic(target, x, cpc)

                if min_pval < self._alpha and (
                    best_var == -1
                    or min_pval < best_min_assoc
                    or (min_pval == best_min_assoc and max_g > best_g_stat)
                ):
                        best_var = x
                        best_min_assoc = min_pval
                        best_g_stat = max_g

            if best_var == -1:
                break

            cpc.append(best_var)
            remaining.remove(best_var)

        return cpc

    def _backward(self, target: int, cpc: list[int]) -> list[int]:
        """Backward phase: remove variables that become independent given subsets of CPC."""
        if len(cpc) <= 1:
            return cpc

        to_remove: list[int] = []
        for x in cpc:
            others = [v for v in cpc if v != x]
            # Test X _||_ T | subsets of CPC \ {X}
            for size in range(len(others) + 1):
                for subset in combinations(others, size):
                    z = list(subset)
                    pval, _ = g_test(self._data, x, target, z, self._cardinalities)
                    if pval > self._alpha and pval != 1.0:
                        to_remove.append(x)
                        break
                else:
                    continue
                break

        return [v for v in cpc if v not in to_remove]

    def _max_min_heuristic(
        self, target: int, x: int, cpc: list[int]
    ) -> tuple[float, float]:
        """Compute the max-min association between x and target given CPC.

        Tests X _||_ T | S for all subsets S of CPC.
        Returns (min_p_value, corresponding_g_statistic) across all subsets.
        The "min" p-value represents the strongest evidence of dependence.
        """
        min_pval = 1.0
        best_g = 0.0

        if len(cpc) == 0:
            pval, g_stat = g_test(self._data, x, target, [], self._cardinalities)
            return pval, g_stat

        # Test with all subsets of CPC as conditioning set
        for size in range(len(cpc) + 1):
            for subset in combinations(cpc, size):
                z = list(subset)
                pval, g_stat = g_test(self._data, x, target, z, self._cardinalities)
                if pval < min_pval or (pval == min_pval and g_stat > best_g):
                    min_pval = pval
                    best_g = g_stat

        return min_pval, best_g

    def _symmetrize(self, pc_sets: list[list[int]]) -> list[list[int]]:
        """Enforce symmetry: x in PC(t) iff t in PC(x)."""
        result: list[list[int]] = [[] for _ in range(self._n_vars)]
        for t in range(self._n_vars):
            for x in pc_sets[t]:
                if t in pc_sets[x]:
                    if x not in result[t]:
                        result[t].append(x)
                    if t not in result[x]:
                        result[x].append(t)
        return result
