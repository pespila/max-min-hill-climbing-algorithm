"""Edge case tests for the MMHC algorithm."""

from __future__ import annotations

import numpy as np

from mmhc import MMHCConfig, mmhc
from mmhc.graph_utils import is_dag


class TestEdgeCases:
    def test_single_variable(self):
        """Single variable should produce empty graph."""
        data = np.array([[0], [1], [0], [1], [0]])
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.shape == (1, 1)
        assert result.adjacency_matrix.sum() == 0

    def test_two_variables_independent(self):
        """Two independent variables should have no edges."""
        rng = np.random.default_rng(42)
        n = 2000
        data = np.column_stack([rng.choice(2, size=n), rng.choice(2, size=n)])
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.sum() == 0

    def test_two_variables_dependent(self):
        """Two dependent variables should have one edge."""
        rng = np.random.default_rng(42)
        n = 2000
        x = rng.choice(2, size=n)
        y = x.copy()
        data = np.column_stack([x, y])
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.sum() == 1
        assert is_dag(result.adjacency_matrix)

    def test_all_identical_values(self):
        """All identical values should produce empty graph."""
        data = np.zeros((100, 3), dtype=np.int64)
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.sum() == 0

    def test_small_sample(self):
        """Small sample should not crash."""
        data = np.array([[0, 1, 0], [1, 0, 1], [0, 0, 1], [1, 1, 0], [0, 1, 1]])
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.shape == (3, 3)
        assert is_dag(result.adjacency_matrix)
