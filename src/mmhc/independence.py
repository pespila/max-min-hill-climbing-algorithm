"""Conditional independence testing via the G-test (chi-squared).

Replaces the C++ Svalue function and D1-D5 multi-dimensional arrays with a
vectorized approach using np.ravel_multi_index + np.bincount, supporting
arbitrary numbers of conditioning variables.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np
from scipy.stats import chi2

if TYPE_CHECKING:
    from numpy.typing import NDArray


@numba.njit(cache=True)
def _g_statistic(obs: NDArray[np.float64], exp: NDArray[np.float64]) -> float:
    """Compute the G-statistic: 2 * sum(obs * ln(obs / exp)) for obs > 0 and exp > 0."""
    g = 0.0
    for i in range(obs.shape[0]):
        if obs[i] > 0.0 and exp[i] > 0.0:
            g += obs[i] * np.log(obs[i] / exp[i])
    return 2.0 * g


def g_test(
    data: NDArray[np.int_],
    x: int,
    y: int,
    z: list[int],
    cardinalities: NDArray[np.int_],
) -> tuple[float, float]:
    """Perform a G-test of conditional independence: X _||_ Y | Z.

    Args:
        data: (n_samples, n_vars) integer matrix, 0-indexed values.
        x: Index of first variable.
        y: Index of second variable.
        z: Indices of conditioning variables (may be empty).
        cardinalities: Number of unique values per variable.

    Returns:
        (p_value, g_statistic) tuple.
    """
    n_samples = data.shape[0]
    rx = int(cardinalities[x])
    ry = int(cardinalities[y])

    if len(z) == 0:
        # Unconditional test: X _||_ Y
        # Build 2D contingency table
        joint_idx = data[:, x] * ry + data[:, y]
        joint_counts = np.bincount(joint_idx, minlength=rx * ry).astype(np.float64)
        joint_table = joint_counts.reshape(rx, ry)

        margin_x = joint_table.sum(axis=1)
        margin_y = joint_table.sum(axis=0)

        expected = np.outer(margin_x, margin_y) / n_samples
        g_stat = _g_statistic(joint_table.ravel(), expected.ravel())

        df = (rx - 1) * (ry - 1)
    else:
        # Conditional test: X _||_ Y | Z
        z_cards = [int(cardinalities[zi]) for zi in z]

        # Encode Z into a single index
        z_cols = np.column_stack([data[:, zi] for zi in z])
        z_strides = np.ones(len(z), dtype=np.int64)
        for i in range(len(z) - 2, -1, -1):
            z_strides[i] = z_strides[i + 1] * z_cards[i + 1]
        z_idx = (z_cols * z_strides).sum(axis=1)

        rz = int(np.prod(np.array(z_cards)))

        # Build 3D contingency table as flat array: (x, y, z)
        flat_idx = (data[:, x] * ry + data[:, y]) * rz + z_idx
        total_cells = rx * ry * rz
        counts = np.bincount(flat_idx, minlength=total_cells).astype(np.float64)
        table_xyz = counts.reshape(rx, ry, rz)

        margin_xz = table_xyz.sum(axis=1)  # (rx, rz)
        margin_yz = table_xyz.sum(axis=0)  # (ry, rz)
        margin_z = table_xyz.sum(axis=(0, 1))  # (rz,)

        # G-statistic: sum over all cells
        g_stat = 0.0
        obs_flat = table_xyz.ravel()
        exp_flat = np.empty_like(obs_flat)
        idx = 0
        for xi in range(rx):
            for yi in range(ry):
                for zi_val in range(rz):
                    mz = margin_z[zi_val]
                    if mz > 0:
                        exp_flat[idx] = margin_xz[xi, zi_val] * margin_yz[yi, zi_val] / mz
                    else:
                        exp_flat[idx] = 0.0
                    idx += 1

        g_stat = _g_statistic(obs_flat, exp_flat)
        df = (rx - 1) * (ry - 1) * rz

    if df <= 0:
        return 1.0, 0.0

    p_value = float(chi2.sf(g_stat, df))
    return p_value, float(g_stat)
