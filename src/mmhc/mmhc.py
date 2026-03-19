"""Top-level MMHC orchestrator class.

Provides a scikit-learn-style fit() API and a convenience mmhc() function.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._types import MMHCConfig, MMHCResult
from .hillclimb import HillClimber
from .mmpc import MMPC

if TYPE_CHECKING:
    from numpy.typing import NDArray


class MMHC:
    """Max-Min Hill-Climbing Bayesian network structure learner.

    Example usage::

        from mmhc import MMHC, make_student

        data = make_student(5000, random_state=42)
        model = MMHC()
        result = model.fit(data)
        print(result.adjacency_matrix)
    """

    def __init__(self, config: MMHCConfig | None = None) -> None:
        self._config = config or MMHCConfig()
        self._result: MMHCResult | None = None

    @property
    def result(self) -> MMHCResult | None:
        """The result of the last fit() call, or None."""
        return self._result

    @property
    def adjacency_matrix_(self) -> NDArray[np.int_]:
        """The learned adjacency matrix. Raises if fit() has not been called."""
        if self._result is None:
            raise ValueError("Model has not been fitted yet. Call fit() first.")
        return self._result.adjacency_matrix

    @property
    def score_(self) -> float:
        """The total BDeu score. Raises if fit() has not been called."""
        if self._result is None:
            raise ValueError("Model has not been fitted yet. Call fit() first.")
        return self._result.score

    def fit(self, data: pd.DataFrame | NDArray[np.int_]) -> MMHCResult:
        """Learn Bayesian network structure from data.

        Args:
            data: DataFrame or integer numpy array. If DataFrame, columns are
                  auto-encoded to 0-indexed integers.

        Returns:
            MMHCResult with the learned DAG.
        """
        int_data, column_names, cardinalities = self._prepare_data(data)

        rng = np.random.default_rng(self._config.random_seed)

        # Phase 1: MMPC skeleton learning
        mmpc = MMPC(int_data, cardinalities, alpha=self._config.alpha)
        pc_sets = mmpc.run()

        # Phase 2: Hill-climbing edge direction
        climber = HillClimber(
            data=int_data,
            cardinalities=cardinalities,
            pc_sets=pc_sets,
            eta=self._config.eta,
            max_iterations=self._config.max_iterations,
            early_stop_rounds=self._config.early_stop_rounds,
            rng=rng,
        )
        adj, node_scores, total_score, n_iters, converged = climber.run()

        self._result = MMHCResult(
            adjacency_matrix=adj,
            parent_children=pc_sets,
            score=total_score,
            node_scores=node_scores,
            n_iterations=n_iters,
            converged=converged,
            column_names=column_names,
        )
        return self._result

    def fit_predict(self, data: pd.DataFrame | NDArray[np.int_]) -> NDArray[np.int_]:
        """Fit and return the adjacency matrix."""
        result = self.fit(data)
        return result.adjacency_matrix

    @staticmethod
    def _prepare_data(
        data: pd.DataFrame | NDArray[np.int_],
    ) -> tuple[NDArray[np.int_], list[str], NDArray[np.int_]]:
        """Convert input data to 0-indexed integer matrix with cardinalities."""
        if isinstance(data, pd.DataFrame):
            column_names = list(data.columns)
            int_data = np.empty((len(data), len(data.columns)), dtype=np.int64)
            for i, col in enumerate(data.columns):
                vals = data[col].values
                if vals.dtype == object or isinstance(vals.dtype, pd.CategoricalDtype):
                    # Encode string/categorical to 0-indexed integers
                    categories = sorted(set(vals))
                    mapping = {v: idx for idx, v in enumerate(categories)}
                    int_data[:, i] = np.array([mapping[v] for v in vals])
                else:
                    # Assume already integer; shift to 0-indexed if needed
                    min_val = int(vals.min())
                    int_data[:, i] = vals.astype(np.int64) - min_val
        else:
            column_names = [str(i) for i in range(data.shape[1])]
            min_vals = data.min(axis=0)
            int_data = (data - min_vals).astype(np.int64)

        cardinalities = np.array(
            [len(np.unique(int_data[:, i])) for i in range(int_data.shape[1])],
            dtype=np.int64,
        )
        return int_data, column_names, cardinalities


def mmhc(
    data: pd.DataFrame | NDArray[np.int_],
    config: MMHCConfig | None = None,
) -> MMHCResult:
    """Convenience function to run the MMHC algorithm.

    Args:
        data: DataFrame or integer numpy array.
        config: Algorithm configuration. Uses defaults if None.

    Returns:
        MMHCResult with the learned DAG.
    """
    model = MMHC(config=config)
    return model.fit(data)
