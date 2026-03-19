"""Shared fixtures for MMHC tests."""

from __future__ import annotations

import numpy as np
import pytest

from mmhc.datasets import make_rainy, make_student


@pytest.fixture
def student_data():
    """Student dataset with 5000 samples, fixed seed."""
    return make_student(5000, random_state=42)


@pytest.fixture
def rainy_data():
    """Rainy dataset with 5000 samples, fixed seed."""
    return make_rainy(5000, random_state=42)


@pytest.fixture
def student_ground_truth():
    """Ground truth adjacency matrix for the student network.

    Edges: difficulty->grade, intelligence->grade, intelligence->SAT, grade->letter
    Columns: difficulty(0), intelligence(1), SAT(2), grade(3), letter(4)
    """
    adj = np.zeros((5, 5), dtype=np.int64)
    adj[0, 3] = 1  # difficulty -> grade
    adj[1, 3] = 1  # intelligence -> grade
    adj[1, 2] = 1  # intelligence -> SAT
    adj[3, 4] = 1  # grade -> letter
    return adj


@pytest.fixture
def rainy_ground_truth():
    """Ground truth adjacency matrix for the rainy network.

    Edges: rain->sprinkler, rain->grassWet, sprinkler->grassWet
    Columns: sprinkler(0), rain(1), grassWet(2)
    """
    adj = np.zeros((3, 3), dtype=np.int64)
    adj[1, 0] = 1  # rain -> sprinkler
    adj[1, 2] = 1  # rain -> grassWet
    adj[0, 2] = 1  # sprinkler -> grassWet
    return adj
