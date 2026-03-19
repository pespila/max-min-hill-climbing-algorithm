"""Dataclasses and type aliases for the MMHC algorithm."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


@dataclass(frozen=True)
class MMHCConfig:
    """Configuration for the MMHC algorithm.

    Attributes:
        alpha: Significance level for conditional independence tests.
        eta: Equivalent sample size for BDeu scoring.
        max_iterations: Maximum hill-climbing iterations.
        early_stop_rounds: Stop after this many non-improving rounds.
        random_seed: Seed for reproducibility. None for non-deterministic.
    """

    alpha: float = 0.05
    eta: float = 1.0
    max_iterations: int = 100
    early_stop_rounds: int = 5
    random_seed: int | None = None


@dataclass
class MMHCResult:
    """Result of the MMHC algorithm.

    Attributes:
        adjacency_matrix: (n_vars, n_vars) binary matrix. Entry [i, j] = 1 means edge i -> j.
        parent_children: List of parent-children sets per variable (0-indexed).
        score: Total BDeu score of the learned DAG.
        node_scores: Per-node BDeu scores.
        n_iterations: Number of hill-climbing iterations performed.
        converged: Whether early stopping was triggered.
        column_names: Original column names from the input DataFrame.
    """

    adjacency_matrix: NDArray[np.int_]
    parent_children: list[list[int]]
    score: float
    node_scores: NDArray[np.float64]
    n_iterations: int
    converged: bool
    column_names: list[str] = field(default_factory=list)
