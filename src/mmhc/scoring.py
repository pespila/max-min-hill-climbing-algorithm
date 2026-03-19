"""BDeu (Bayesian Dirichlet likelihood-equivalence uniform) scoring.

Replaces the C++ ScoreNodeWith{None,One,More}Parents functions with a unified
vectorized implementation.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numba
import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


@numba.njit(cache=True)
def _bdeu_score_counts(
    n_ij: NDArray[np.float64],
    n_ijk: NDArray[np.float64],
    q: int,
    r: int,
    eta: float,
) -> float:
    """Compute BDeu score from count arrays.

    Args:
        n_ij: (q,) array of parent-config counts.
        n_ijk: (q, r) array of joint counts.
        q: Number of parent configurations.
        r: Cardinality of the child node.
        eta: Equivalent sample size.

    Returns:
        BDeu score (log-scale).
    """
    gamma_j = eta / q
    gamma_k = eta / (q * r)
    score = 0.0

    for j in range(q):
        score += math.lgamma(gamma_j) - math.lgamma(n_ij[j] + gamma_j)
        for k in range(r):
            if n_ijk[j, k] > 0:
                score += math.lgamma(n_ijk[j, k] + gamma_k) - math.lgamma(gamma_k)

    return score


def bdeu_score_node(
    data: NDArray[np.int_],
    node: int,
    parents: list[int],
    cardinalities: NDArray[np.int_],
    eta: float = 1.0,
) -> float:
    """Compute the BDeu score for a single node given its parents.

    Args:
        data: (n_samples, n_vars) integer matrix, 0-indexed values.
        node: Index of the child node.
        parents: Indices of parent nodes (may be empty).
        cardinalities: Number of unique values per variable.
        eta: Equivalent sample size.

    Returns:
        BDeu score (log-scale, higher is better).
    """
    r = int(cardinalities[node])
    n_samples = data.shape[0]

    if len(parents) == 0:
        q = 1
        child_vals = data[:, node]
        n_ijk = np.bincount(child_vals, minlength=r).astype(np.float64).reshape(1, r)
        n_ij = np.array([float(n_samples)])
    else:
        parent_cards = [int(cardinalities[p]) for p in parents]
        q = int(np.prod(np.array(parent_cards)))

        # Encode parent configurations into a single index
        parent_cols = data[:, parents]
        strides = np.ones(len(parents), dtype=np.int64)
        for i in range(len(parents) - 2, -1, -1):
            strides[i] = strides[i + 1] * parent_cards[i + 1]
        parent_idx = (parent_cols * strides).sum(axis=1)

        # Build joint count table (q, r)
        flat_idx = parent_idx * r + data[:, node]
        counts = np.bincount(flat_idx, minlength=q * r).astype(np.float64)
        n_ijk = counts.reshape(q, r)
        n_ij = n_ijk.sum(axis=1)

    return float(_bdeu_score_counts(n_ij, n_ijk, q, r, eta))


def bdeu_score_graph(
    data: NDArray[np.int_],
    adjacency: NDArray[np.int_],
    cardinalities: NDArray[np.int_],
    eta: float = 1.0,
) -> NDArray[np.float64]:
    """Compute per-node BDeu scores for an entire graph.

    Args:
        data: (n_samples, n_vars) integer matrix, 0-indexed values.
        adjacency: (n_vars, n_vars) binary matrix. Entry [i, j] = 1 means i -> j.
        cardinalities: Number of unique values per variable.
        eta: Equivalent sample size.

    Returns:
        (n_vars,) array of per-node BDeu scores.
    """
    n_vars = data.shape[1]
    scores = np.empty(n_vars, dtype=np.float64)
    for node in range(n_vars):
        parents = list(np.where(adjacency[:, node] == 1)[0])
        scores[node] = bdeu_score_node(data, node, parents, cardinalities, eta)
    return scores
