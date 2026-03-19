"""Tests for the hill-climbing edge direction module."""

from __future__ import annotations

import numpy as np

from mmhc.graph_utils import is_dag
from mmhc.hillclimb import HillClimber


class TestHillClimber:
    def test_produces_dag(self):
        """Hill-climbing should always produce a DAG (no cycles)."""
        rng = np.random.default_rng(42)
        n = 1000
        x = rng.choice(2, size=n)
        y = x.copy()
        z = rng.choice(2, size=n)
        data = np.column_stack([x, y, z])
        cardinalities = np.array([2, 2, 2])

        # Give it a PC skeleton
        pc_sets = [[1], [0], []]
        climber = HillClimber(
            data, cardinalities, pc_sets, eta=1.0,
            max_iterations=20, early_stop_rounds=3, rng=np.random.default_rng(42),
        )
        adj, scores, total, n_iters, converged = climber.run()
        assert is_dag(adj)

    def test_returns_correct_shape(self):
        """Output shapes should match number of variables."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(500, 4))
        cardinalities = np.array([2, 2, 2, 2])
        pc_sets = [[1], [0, 2], [1, 3], [2]]

        climber = HillClimber(
            data, cardinalities, pc_sets, eta=1.0,
            max_iterations=10, early_stop_rounds=3, rng=np.random.default_rng(42),
        )
        adj, scores, total, n_iters, converged = climber.run()
        assert adj.shape == (4, 4)
        assert scores.shape == (4,)
        assert isinstance(total, float)
        assert isinstance(n_iters, int)
        assert isinstance(converged, bool)

    def test_deterministic_with_seed(self):
        """Same seed should produce identical results."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(500, 3))
        cardinalities = np.array([2, 2, 2])
        pc_sets = [[1], [0, 2], [1]]

        results = []
        for _ in range(2):
            climber = HillClimber(
                data, cardinalities, pc_sets, eta=1.0,
                max_iterations=10, early_stop_rounds=3, rng=np.random.default_rng(99),
            )
            adj, _, total, _, _ = climber.run()
            results.append((adj.copy(), total))

        np.testing.assert_array_equal(results[0][0], results[1][0])
        assert results[0][1] == results[1][1]

    def test_empty_skeleton(self):
        """Empty PC sets should produce empty graph."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(200, 3))
        cardinalities = np.array([2, 2, 2])
        pc_sets: list[list[int]] = [[], [], []]

        climber = HillClimber(
            data, cardinalities, pc_sets, eta=1.0,
            max_iterations=5, early_stop_rounds=3, rng=np.random.default_rng(42),
        )
        adj, _, _, _, _ = climber.run()
        assert adj.sum() == 0
