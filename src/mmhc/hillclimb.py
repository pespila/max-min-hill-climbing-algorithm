"""Hill-climbing edge direction using BDeu scoring.

Replaces the C++ SettingEdges and AddReverseDelete functions.
Adds cycle checking via graph_utils.would_create_cycle().
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .graph_utils import would_create_cycle
from .scoring import bdeu_score_graph

if TYPE_CHECKING:
    from numpy.typing import NDArray


class HillClimber:
    """Greedy hill-climbing for directing edges in a Bayesian network skeleton."""

    def __init__(
        self,
        data: NDArray[np.int_],
        cardinalities: NDArray[np.int_],
        pc_sets: list[list[int]],
        eta: float = 1.0,
        max_iterations: int = 100,
        early_stop_rounds: int = 5,
        rng: np.random.Generator | None = None,
    ) -> None:
        self._data = data
        self._cardinalities = cardinalities
        self._pc_sets = pc_sets
        self._eta = eta
        self._max_iterations = max_iterations
        self._early_stop_rounds = early_stop_rounds
        self._rng = rng if rng is not None else np.random.default_rng()
        self._n_vars = data.shape[1]

    def run(
        self,
    ) -> tuple[NDArray[np.int_], NDArray[np.float64], float, int, bool]:
        """Run hill-climbing to direct edges.

        Returns:
            (adjacency_matrix, node_scores, total_score, n_iterations, converged)
        """
        adj = np.zeros((self._n_vars, self._n_vars), dtype=np.int_)
        scores = bdeu_score_graph(self._data, adj, self._cardinalities, self._eta)

        # Phase 1: Set initial edges from skeleton
        self._set_initial_edges(adj, scores)

        # Phase 2: Iterative add/reverse/delete
        no_improve_count = 0
        n_iters = 0
        converged = False

        for i in range(self._max_iterations):
            n_iters = i + 1
            old_total = float(scores.sum())

            tmp_adj = adj.copy()
            tmp_scores = scores.copy()
            self._add_reverse_delete(tmp_adj, tmp_scores)

            new_total = float(tmp_scores.sum())
            if new_total > old_total:
                adj[:] = tmp_adj
                scores[:] = tmp_scores
                no_improve_count = 0
            else:
                no_improve_count += 1

            if no_improve_count >= self._early_stop_rounds:
                converged = True
                break

        total_score = float(scores.sum())
        return adj, scores, total_score, n_iters, converged

    def _set_initial_edges(
        self, adj: NDArray[np.int_], scores: NDArray[np.float64]
    ) -> None:
        """Initialize edges from the PC skeleton, keeping only score-improving ones."""
        for i in range(self._n_vars):
            for j in self._pc_sets[i]:
                if adj[i, j] == 0 and adj[j, i] == 0:
                    # Try edge i -> j
                    if would_create_cycle(adj, i, j):
                        continue
                    adj[i, j] = 1
                    new_scores = bdeu_score_graph(
                        self._data, adj, self._cardinalities, self._eta
                    )
                    if float(new_scores.sum()) > float(scores.sum()):
                        scores[:] = new_scores
                    else:
                        adj[i, j] = 0

    def _add_reverse_delete(
        self, adj: NDArray[np.int_], scores: NDArray[np.float64]
    ) -> None:
        """Try add, reverse, and delete operations on edges."""
        for i in range(self._n_vars):
            for j in self._pc_sets[i]:
                rnd = self._rng.integers(1, 11)

                if adj[i, j] == 0 and adj[j, i] == 0:
                    # Add: try reverse direction j -> i
                    if would_create_cycle(adj, j, i):
                        continue
                    adj[j, i] = 1
                    self._accept_or_revert(adj, scores, j, i)

                elif adj[i, j] == 1 and rnd > 5:
                    # Delete edge i -> j
                    adj[i, j] = 0
                    self._accept_or_revert_delete(adj, scores, i, j)

                elif adj[j, i] == 1 and rnd <= 5:
                    # Reverse edge j -> i to i -> j
                    adj[j, i] = 0
                    if would_create_cycle(adj, i, j):
                        adj[j, i] = 1
                        continue
                    adj[i, j] = 1
                    new_scores = bdeu_score_graph(
                        self._data, adj, self._cardinalities, self._eta
                    )
                    if float(new_scores.sum()) > float(scores.sum()):
                        scores[:] = new_scores
                    else:
                        adj[i, j] = 0
                        adj[j, i] = 1

    def _accept_or_revert(
        self,
        adj: NDArray[np.int_],
        scores: NDArray[np.float64],
        src: int,
        dst: int,
    ) -> None:
        """Accept edge addition if score improves, otherwise revert."""
        new_scores = bdeu_score_graph(
            self._data, adj, self._cardinalities, self._eta
        )
        if float(new_scores.sum()) > float(scores.sum()):
            scores[:] = new_scores
        else:
            adj[src, dst] = 0

    def _accept_or_revert_delete(
        self,
        adj: NDArray[np.int_],
        scores: NDArray[np.float64],
        src: int,
        dst: int,
    ) -> None:
        """Accept edge deletion if score improves, otherwise revert."""
        new_scores = bdeu_score_graph(
            self._data, adj, self._cardinalities, self._eta
        )
        if float(new_scores.sum()) > float(scores.sum()):
            scores[:] = new_scores
        else:
            adj[src, dst] = 1
