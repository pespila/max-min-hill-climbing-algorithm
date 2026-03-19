"""MMHC — Max-Min Hill-Climbing Bayesian network structure learning."""

from ._types import MMHCConfig, MMHCResult
from .datasets import make_rainy, make_student
from .mmhc import MMHC, mmhc
from .plotting import plot_dag

__version__ = "2.0.0"

__all__ = [
    "MMHC",
    "MMHCConfig",
    "MMHCResult",
    "make_rainy",
    "make_student",
    "mmhc",
    "plot_dag",
]
