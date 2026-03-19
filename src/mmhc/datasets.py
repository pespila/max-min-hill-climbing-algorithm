"""Synthetic dataset generators for testing the MMHC algorithm.

Reproduces the exact CPTs from the original R implementation in data.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def make_student(n_samples: int = 1000, random_state: int | None = None) -> pd.DataFrame:
    """Generate synthetic student data from a known Bayesian network.

    Network structure:
        difficulty -> grade
        intelligence -> grade
        intelligence -> SAT
        grade -> letter

    All variables are 0-indexed integers:
        difficulty: {0, 1}
        intelligence: {0, 1}
        SAT: {0, 1}
        grade: {0, 1, 2}
        letter: {0, 1}

    The CPTs match the original R implementation (shifted from 1-indexed to 0-indexed).
    """
    rng = np.random.default_rng(random_state)

    difficulty = rng.choice(2, size=n_samples, p=[0.6, 0.4])
    intelligence = rng.choice(2, size=n_samples, p=[0.7, 0.3])

    sat = np.empty(n_samples, dtype=np.int64)
    grade = np.empty(n_samples, dtype=np.int64)
    letter = np.empty(n_samples, dtype=np.int64)

    for n in range(n_samples):
        i = intelligence[n]
        d = difficulty[n]

        # SAT depends on intelligence
        if i == 0:
            sat[n] = rng.choice(2, p=[0.95, 0.05])
        else:
            sat[n] = rng.choice(2, p=[0.2, 0.8])

        # grade depends on intelligence and difficulty
        if i == 1 and d == 1:
            grade[n] = rng.choice(3, p=[0.5, 0.3, 0.2])
        elif i == 1 and d == 0:
            grade[n] = rng.choice(3, p=[0.9, 0.08, 0.02])
        elif i == 0 and d == 1:
            grade[n] = rng.choice(3, p=[0.05, 0.25, 0.7])
        else:  # i == 0 and d == 0
            grade[n] = rng.choice(3, p=[0.3, 0.4, 0.3])

        # letter depends on grade
        g = grade[n]
        if g == 0:
            letter[n] = rng.choice(2, p=[0.1, 0.9])
        elif g == 1:
            letter[n] = rng.choice(2, p=[0.4, 0.6])
        else:  # g == 2
            letter[n] = rng.choice(2, p=[0.99, 0.01])

    return pd.DataFrame({
        "difficulty": difficulty,
        "intelligence": intelligence,
        "SAT": sat,
        "grade": grade,
        "letter": letter,
    })


def make_rainy(n_samples: int = 1000, random_state: int | None = None) -> pd.DataFrame:
    """Generate synthetic weather/sprinkler data from a known Bayesian network.

    Network structure:
        rain -> sprinkler
        rain -> grassWet
        sprinkler -> grassWet

    All variables are 0-indexed integers: {0, 1}.

    The CPTs match the original R implementation (shifted from 1-indexed to 0-indexed).
    """
    rng = np.random.default_rng(random_state)

    rain = rng.choice(2, size=n_samples, p=[0.2, 0.8])

    sprinkler = np.empty(n_samples, dtype=np.int64)
    grass_wet = np.empty(n_samples, dtype=np.int64)

    for n in range(n_samples):
        r = rain[n]

        # sprinkler depends on rain
        if r == 0:
            sprinkler[n] = rng.choice(2, p=[0.01, 0.99])
        else:
            sprinkler[n] = rng.choice(2, p=[0.40, 0.60])

        s = sprinkler[n]

        # grassWet depends on rain and sprinkler
        if r == 0 and s == 0:
            grass_wet[n] = rng.choice(2, p=[0.99, 0.01])
        elif r == 0 and s == 1:
            grass_wet[n] = rng.choice(2, p=[0.80, 0.20])
        elif r == 1 and s == 0:
            grass_wet[n] = rng.choice(2, p=[0.90, 0.10])
        else:  # r == 1 and s == 1
            grass_wet[n] = rng.choice(2, p=[0.00, 1.00])

    return pd.DataFrame({
        "sprinkler": sprinkler,
        "rain": rain,
        "grassWet": grass_wet,
    })
