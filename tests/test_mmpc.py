"""Tests for the MMPC skeleton learning module."""

from __future__ import annotations

import numpy as np

from mmhc.datasets import make_rainy, make_student
from mmhc.mmpc import MMPC


class TestMMPC:
    def test_student_skeleton(self):
        """MMPC should recover the student network skeleton."""
        df = make_student(5000, random_state=42)
        data = df.values
        min_vals = data.min(axis=0)
        int_data = (data - min_vals).astype(np.int64)
        cardinalities = np.array([len(np.unique(int_data[:, i])) for i in range(int_data.shape[1])])

        mmpc = MMPC(int_data, cardinalities, alpha=0.05)
        pc_sets = mmpc.run()

        assert len(pc_sets) == 5
        # difficulty(0) should be connected to grade(3)
        assert 3 in pc_sets[0]
        # intelligence(1) should be connected to grade(3) and SAT(2)
        assert 3 in pc_sets[1]
        assert 2 in pc_sets[1]
        # grade(3) should be connected to letter(4)
        assert 4 in pc_sets[3]

    def test_rainy_skeleton(self):
        """MMPC should recover the rainy network skeleton."""
        df = make_rainy(5000, random_state=42)
        data = df.values
        min_vals = data.min(axis=0)
        int_data = (data - min_vals).astype(np.int64)
        cardinalities = np.array([len(np.unique(int_data[:, i])) for i in range(int_data.shape[1])])

        mmpc = MMPC(int_data, cardinalities, alpha=0.05)
        pc_sets = mmpc.run()

        assert len(pc_sets) == 3
        # rain(1) should be connected to sprinkler(0) and grassWet(2)
        assert 0 in pc_sets[1]
        assert 2 in pc_sets[1]

    def test_symmetry(self):
        """PC sets should be symmetric: x in PC(t) iff t in PC(x)."""
        df = make_student(3000, random_state=42)
        data = df.values
        min_vals = data.min(axis=0)
        int_data = (data - min_vals).astype(np.int64)
        cardinalities = np.array([len(np.unique(int_data[:, i])) for i in range(int_data.shape[1])])

        mmpc = MMPC(int_data, cardinalities, alpha=0.05)
        pc_sets = mmpc.run()

        for t, pc in enumerate(pc_sets):
            for x in pc:
                assert t in pc_sets[x], f"{t} in PC({x}) should hold since {x} in PC({t})"

    def test_independent_variables(self):
        """Independent variables should have empty PC sets."""
        rng = np.random.default_rng(42)
        n = 2000
        data = np.column_stack([
            rng.choice(2, size=n),
            rng.choice(3, size=n),
            rng.choice(2, size=n),
        ])
        cardinalities = np.array([2, 3, 2])

        mmpc = MMPC(data, cardinalities, alpha=0.05)
        pc_sets = mmpc.run()

        for pc in pc_sets:
            assert len(pc) == 0
