"""Tests for the BDeu scoring module."""

from __future__ import annotations

import numpy as np

from mmhc.scoring import bdeu_score_graph, bdeu_score_node


class TestBdeuScoreNode:
    def test_no_parents(self):
        """Score with no parents should be computable."""
        rng = np.random.default_rng(42)
        data = rng.choice(3, size=(100, 3))
        cardinalities = np.array([3, 3, 3])
        score = bdeu_score_node(data, 0, [], cardinalities)
        assert isinstance(score, float)
        assert score < 0  # Log-likelihood is always negative

    def test_one_parent(self):
        """Score with one parent."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(500, 3))
        cardinalities = np.array([2, 2, 2])
        score = bdeu_score_node(data, 1, [0], cardinalities)
        assert isinstance(score, float)

    def test_multiple_parents(self):
        """Score with multiple parents."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(500, 4))
        cardinalities = np.array([2, 2, 2, 2])
        score = bdeu_score_node(data, 3, [0, 1], cardinalities)
        assert isinstance(score, float)

    def test_true_parent_improves_score(self):
        """Adding a true parent should improve the score."""
        rng = np.random.default_rng(42)
        n = 2000
        x = rng.choice(2, size=n)
        y = x.copy()  # y is deterministic function of x
        data = np.column_stack([x, y])
        cardinalities = np.array([2, 2])

        score_no_parent = bdeu_score_node(data, 1, [], cardinalities)
        score_with_parent = bdeu_score_node(data, 1, [0], cardinalities)
        assert score_with_parent > score_no_parent

    def test_false_parent_worsens_score(self):
        """Adding an independent variable as parent should worsen score."""
        rng = np.random.default_rng(42)
        n = 2000
        data = rng.choice(2, size=(n, 2))  # Independent variables
        cardinalities = np.array([2, 2])

        score_no_parent = bdeu_score_node(data, 1, [], cardinalities)
        score_with_parent = bdeu_score_node(data, 1, [0], cardinalities)
        assert score_with_parent < score_no_parent


class TestBdeuScoreGraph:
    def test_empty_graph(self):
        """Score of empty graph should be sum of no-parent scores."""
        rng = np.random.default_rng(42)
        data = rng.choice(2, size=(200, 3))
        adj = np.zeros((3, 3), dtype=np.int64)
        cardinalities = np.array([2, 2, 2])

        scores = bdeu_score_graph(data, adj, cardinalities)
        assert scores.shape == (3,)
        for i in range(3):
            expected = bdeu_score_node(data, i, [], cardinalities)
            assert abs(scores[i] - expected) < 1e-10

    def test_graph_with_edges(self):
        """Score with edges should differ from empty graph."""
        rng = np.random.default_rng(42)
        n = 500
        x = rng.choice(2, size=n)
        y = x.copy()
        data = np.column_stack([x, y])
        cardinalities = np.array([2, 2])

        adj_empty = np.zeros((2, 2), dtype=np.int64)
        adj_edge = np.zeros((2, 2), dtype=np.int64)
        adj_edge[0, 1] = 1  # x -> y

        scores_empty = bdeu_score_graph(data, adj_empty, cardinalities)
        scores_edge = bdeu_score_graph(data, adj_edge, cardinalities)
        assert scores_edge.sum() > scores_empty.sum()
