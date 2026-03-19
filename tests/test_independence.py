"""Tests for the G-test conditional independence module."""

from __future__ import annotations

import numpy as np

from mmhc.independence import g_test


class TestGTestUnconditional:
    def test_independent_variables(self):
        """Two independent variables should have high p-value."""
        rng = np.random.default_rng(42)
        n = 5000
        data = np.column_stack([
            rng.choice(2, size=n),
            rng.choice(3, size=n),
        ])
        cardinalities = np.array([2, 3])
        pval, g_stat = g_test(data, 0, 1, [], cardinalities)
        assert pval > 0.05

    def test_dependent_variables(self):
        """Two perfectly dependent variables should have low p-value."""
        rng = np.random.default_rng(42)
        n = 1000
        x = rng.choice(3, size=n)
        data = np.column_stack([x, x])
        cardinalities = np.array([3, 3])
        pval, g_stat = g_test(data, 0, 1, [], cardinalities)
        assert pval < 0.001
        assert g_stat > 0

    def test_returns_tuple(self):
        data = np.array([[0, 0], [1, 1], [0, 1], [1, 0]])
        cardinalities = np.array([2, 2])
        result = g_test(data, 0, 1, [], cardinalities)
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestGTestConditional:
    def test_conditional_independence(self):
        """X _||_ Y | Z when X and Y are both caused by Z."""
        rng = np.random.default_rng(42)
        n = 5000
        z = rng.choice(2, size=n)
        x = np.where(z == 0, rng.choice(2, size=n, p=[0.9, 0.1]),
                     rng.choice(2, size=n, p=[0.1, 0.9]))
        y = np.where(z == 0, rng.choice(2, size=n, p=[0.8, 0.2]),
                     rng.choice(2, size=n, p=[0.2, 0.8]))
        data = np.column_stack([x, y, z])
        cardinalities = np.array([2, 2, 2])
        pval, _ = g_test(data, 0, 1, [2], cardinalities)
        assert pval > 0.05

    def test_conditional_dependence(self):
        """X and Y remain dependent even conditioning on Z."""
        rng = np.random.default_rng(42)
        n = 5000
        x = rng.choice(2, size=n)
        y = x.copy()  # Perfect dependence
        z = rng.choice(2, size=n)
        data = np.column_stack([x, y, z])
        cardinalities = np.array([2, 2, 2])
        pval, _ = g_test(data, 0, 1, [2], cardinalities)
        assert pval < 0.001

    def test_multiple_conditioning_vars(self):
        """Test with 2 conditioning variables."""
        rng = np.random.default_rng(42)
        n = 5000
        data = np.column_stack([
            rng.choice(2, size=n),
            rng.choice(2, size=n),
            rng.choice(2, size=n),
            rng.choice(2, size=n),
        ])
        cardinalities = np.array([2, 2, 2, 2])
        pval, _ = g_test(data, 0, 1, [2, 3], cardinalities)
        assert pval > 0.01  # Should be independent
