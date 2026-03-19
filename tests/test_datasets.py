"""Tests for dataset generators."""

from __future__ import annotations

from mmhc.datasets import make_rainy, make_student


class TestMakeStudent:
    def test_shape(self):
        df = make_student(100, random_state=0)
        assert df.shape == (100, 5)

    def test_columns(self):
        df = make_student(10, random_state=0)
        assert list(df.columns) == ["difficulty", "intelligence", "SAT", "grade", "letter"]

    def test_value_ranges(self):
        df = make_student(1000, random_state=0)
        assert set(df["difficulty"].unique()).issubset({0, 1})
        assert set(df["intelligence"].unique()).issubset({0, 1})
        assert set(df["SAT"].unique()).issubset({0, 1})
        assert set(df["grade"].unique()).issubset({0, 1, 2})
        assert set(df["letter"].unique()).issubset({0, 1})

    def test_reproducibility(self):
        df1 = make_student(100, random_state=42)
        df2 = make_student(100, random_state=42)
        assert df1.equals(df2)

    def test_different_seeds(self):
        df1 = make_student(100, random_state=1)
        df2 = make_student(100, random_state=2)
        assert not df1.equals(df2)

    def test_approximate_marginals(self):
        df = make_student(50000, random_state=0)
        # difficulty: P(0) = 0.6
        assert abs(df["difficulty"].mean() - 0.4) < 0.02
        # intelligence: P(0) = 0.7
        assert abs(df["intelligence"].mean() - 0.3) < 0.02


class TestMakeRainy:
    def test_shape(self):
        df = make_rainy(100, random_state=0)
        assert df.shape == (100, 3)

    def test_columns(self):
        df = make_rainy(10, random_state=0)
        assert list(df.columns) == ["sprinkler", "rain", "grassWet"]

    def test_value_ranges(self):
        df = make_rainy(1000, random_state=0)
        for col in df.columns:
            assert set(df[col].unique()).issubset({0, 1})

    def test_reproducibility(self):
        df1 = make_rainy(100, random_state=42)
        df2 = make_rainy(100, random_state=42)
        assert df1.equals(df2)

    def test_approximate_marginals(self):
        df = make_rainy(50000, random_state=0)
        # rain: P(0) = 0.2
        assert abs(df["rain"].mean() - 0.8) < 0.02
