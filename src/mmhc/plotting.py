"""DAG visualization using NetworkX and matplotlib."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .graph_utils import adjacency_to_networkx

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


def plot_dag(
    adjacency: NDArray[np.int_],
    labels: list[str] | None = None,
    ax: Any = None,
    title: str | None = None,
) -> Any:
    """Plot a DAG from an adjacency matrix.

    Args:
        adjacency: (n_vars, n_vars) binary adjacency matrix.
        labels: Node labels. Defaults to "0", "1", ...
        ax: Matplotlib axes. If None, creates a new figure.
        title: Optional plot title.

    Returns:
        The matplotlib axes object.
    """
    import matplotlib.pyplot as plt
    import networkx as nx

    g = adjacency_to_networkx(adjacency, labels)

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 6))

    pos = nx.spring_layout(g, seed=42)
    nx.draw(
        g,
        pos,
        ax=ax,
        with_labels=True,
        node_color="lightblue",
        node_size=2000,
        font_size=10,
        font_weight="bold",
        arrows=True,
        arrowsize=20,
        edge_color="gray",
        width=2,
    )

    if title:
        ax.set_title(title, fontsize=14)

    return ax
