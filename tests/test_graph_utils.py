"""Tests for graph_utils module."""

from __future__ import annotations

import numpy as np

from mmhc.graph_utils import (
    adjacency_to_edges,
    adjacency_to_networkx,
    get_parents,
    is_dag,
    would_create_cycle,
)


class TestGetParents:
    def test_no_parents(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        assert get_parents(0, adj) == []

    def test_single_parent(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1  # 0 -> 1
        assert get_parents(1, adj) == [0]

    def test_multiple_parents(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 2] = 1  # 0 -> 2
        adj[1, 2] = 1  # 1 -> 2
        assert get_parents(2, adj) == [0, 1]


class TestWouldCreateCycle:
    def test_no_cycle_empty(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        assert not would_create_cycle(adj, 0, 1)

    def test_cycle_detected(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1  # 0 -> 1
        adj[1, 2] = 1  # 1 -> 2
        # Adding 2 -> 0 would create cycle
        assert would_create_cycle(adj, 2, 0)

    def test_no_cycle_allowed(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1  # 0 -> 1
        # Adding 0 -> 2 is fine
        assert not would_create_cycle(adj, 0, 2)

    def test_self_loop(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        assert would_create_cycle(adj, 0, 0)


class TestIsDag:
    def test_empty_graph(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        assert is_dag(adj)

    def test_valid_dag(self):
        adj = np.zeros((4, 4), dtype=np.int64)
        adj[0, 1] = 1
        adj[0, 2] = 1
        adj[1, 3] = 1
        adj[2, 3] = 1
        assert is_dag(adj)

    def test_cycle(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1
        adj[1, 2] = 1
        adj[2, 0] = 1
        assert not is_dag(adj)

    def test_single_node(self):
        adj = np.zeros((1, 1), dtype=np.int64)
        assert is_dag(adj)


class TestAdjacencyToEdges:
    def test_empty(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        assert adjacency_to_edges(adj) == []

    def test_edges(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1
        adj[1, 2] = 1
        edges = adjacency_to_edges(adj)
        assert (0, 1) in edges
        assert (1, 2) in edges
        assert len(edges) == 2


class TestAdjacencyToNetworkx:
    def test_basic(self):
        adj = np.zeros((3, 3), dtype=np.int64)
        adj[0, 1] = 1
        g = adjacency_to_networkx(adj, ["A", "B", "C"])
        assert g.has_edge("A", "B")
        assert not g.has_edge("B", "A")
        assert len(g.edges) == 1

    def test_default_labels(self):
        adj = np.zeros((2, 2), dtype=np.int64)
        adj[0, 1] = 1
        g = adjacency_to_networkx(adj)
        assert g.has_edge("0", "1")
