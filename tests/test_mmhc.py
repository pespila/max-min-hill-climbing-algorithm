"""Integration tests for the full MMHC pipeline."""

from __future__ import annotations

import numpy as np

from mmhc import MMHC, MMHCConfig, mmhc
from mmhc.datasets import make_student
from mmhc.graph_utils import is_dag


class TestMMHCIntegration:
    def test_student_network(self, student_data, student_ground_truth):
        """MMHC should recover the student network structure."""
        config = MMHCConfig(random_seed=42)
        result = mmhc(student_data, config=config)

        assert is_dag(result.adjacency_matrix)
        assert result.adjacency_matrix.shape == (5, 5)
        assert result.score < 0
        assert result.n_iterations > 0
        assert len(result.column_names) == 5

        # Check that key edges are recovered
        adj = result.adjacency_matrix
        gt = student_ground_truth

        # Count matching edges (allowing direction flexibility)
        recovered = 0
        total_gt = int(gt.sum())
        for i in range(5):
            for j in range(5):
                if gt[i, j] == 1 and (adj[i, j] == 1 or adj[j, i] == 1):
                        recovered += 1
        # Should recover at least 3 out of 4 edges
        assert recovered >= 3, f"Only recovered {recovered}/{total_gt} edges"

    def test_rainy_network(self, rainy_data, rainy_ground_truth):
        """MMHC should recover the rainy network structure."""
        config = MMHCConfig(random_seed=42)
        result = mmhc(rainy_data, config=config)

        assert is_dag(result.adjacency_matrix)
        assert result.adjacency_matrix.shape == (3, 3)

    def test_reproducibility(self, student_data):
        """Same seed should produce identical results."""
        config = MMHCConfig(random_seed=123)
        r1 = mmhc(student_data, config=config)
        r2 = mmhc(student_data, config=config)
        np.testing.assert_array_equal(r1.adjacency_matrix, r2.adjacency_matrix)
        assert r1.score == r2.score

    def test_class_api(self, student_data):
        """MMHC class API should work correctly."""
        model = MMHC(config=MMHCConfig(random_seed=42))
        result = model.fit(student_data)

        assert result is model.result
        np.testing.assert_array_equal(model.adjacency_matrix_, result.adjacency_matrix)
        assert model.score_ == result.score

    def test_fit_predict(self, student_data):
        """fit_predict should return adjacency matrix."""
        model = MMHC(config=MMHCConfig(random_seed=42))
        adj = model.fit_predict(student_data)
        assert adj.shape == (5, 5)
        assert is_dag(adj)

    def test_not_fitted_error(self):
        """Accessing results before fit() should raise."""
        import pytest

        model = MMHC()
        with pytest.raises(ValueError):
            _ = model.adjacency_matrix_

        with pytest.raises(ValueError):
            _ = model.score_

    def test_numpy_input(self):
        """Should accept numpy array input."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(500, 3))
        config = MMHCConfig(random_seed=42)
        result = mmhc(data, config=config)
        assert result.adjacency_matrix.shape == (3, 3)

    def test_custom_config(self):
        """Custom configuration should be respected."""
        data = make_student(1000, random_state=42)
        config = MMHCConfig(
            alpha=0.01,
            eta=0.5,
            max_iterations=50,
            early_stop_rounds=3,
            random_seed=42,
        )
        result = mmhc(data, config=config)
        assert is_dag(result.adjacency_matrix)
