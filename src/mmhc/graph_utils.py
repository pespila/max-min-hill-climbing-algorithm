"""DAG validation, cycle detection, and adjacency helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


def get_parents(node: int, adjacency: NDArray[np.int_]) -> list[int]:
    """Return the list of parent indices for a given node.

    In the adjacency matrix, entry [i, j] = 1 means edge i -> j,
    so parents of `node` are all `i` where adjacency[i, node] == 1.
    """
    return list(np.where(adjacency[:, node] == 1)[0])


def would_create_cycle(adjacency: NDArray[np.int_], src: int, dst: int) -> bool:
    """Check if adding edge src -> dst would create a cycle.

    Uses DFS from dst to see if src is reachable (which would mean a cycle).
    """
    n = adjacency.shape[0]
    visited = np.zeros(n, dtype=np.bool_)
    stack = [dst]
    while stack:
        current = stack.pop()
        if current == src:
            return True
        if visited[current]:
            continue
        visited[current] = True
        children = np.where(adjacency[current, :] == 1)[0]
        for child in children:
            if not visited[child]:
                stack.append(int(child))
    return False


def is_dag(adjacency: NDArray[np.int_]) -> bool:
    """Check if the adjacency matrix represents a valid DAG (no cycles).

    Uses Kahn's algorithm (topological sort via in-degree counting).
    """
    n = adjacency.shape[0]
    in_degree = adjacency.sum(axis=0).copy()
    queue = list(np.where(in_degree == 0)[0])
    count = 0
    while queue:
        node = queue.pop(0)
        count += 1
        children = np.where(adjacency[node, :] == 1)[0]
        for child in children:
            in_degree[child] -= 1
            if in_degree[child] == 0:
                queue.append(int(child))
    return count == n


def adjacency_to_edges(adjacency: NDArray[np.int_]) -> list[tuple[int, int]]:
    """Convert adjacency matrix to list of (src, dst) edge tuples."""
    rows, cols = np.where(adjacency == 1)
    return list(zip(rows.tolist(), cols.tolist(), strict=True))


def adjacency_to_networkx(
    adjacency: NDArray[np.int_],
    labels: list[str] | None = None,
) -> nx.DiGraph:  # type: ignore[name-defined]  # noqa: F821
    """Convert adjacency matrix to a NetworkX DiGraph."""
    import networkx as nx

    n = adjacency.shape[0]
    if labels is None:
        labels = [str(i) for i in range(n)]

    g = nx.DiGraph()
    g.add_nodes_from(labels)
    for src, dst in adjacency_to_edges(adjacency):
        g.add_edge(labels[src], labels[dst])
    return g
