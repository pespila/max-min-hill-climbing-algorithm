# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python implementation of the **Max-Min Hill-Climbing (MMHC)** Bayesian network structure learning algorithm. Based on the paper by Tsamardinos, Brown & Aliferis (Mach Learn, DOI 10.1007/s10994-006-6889-7).

The algorithm reconstructs Bayesian Networks from observational data in two phases:
1. **MMPC (Max-Min Parents and Children)**: Builds the skeleton (undirected graph) using conditional independence tests (G-test)
2. **MMHC**: Directs edges using greedy hill-climbing with BDeu (Bayesian Dirichlet likelihood-equivalence uniform) scoring

## Build & Install

```bash
cd mmhc
pip install -e ".[dev]"
```

**Dependencies**: Python >= 3.10, numpy, scipy, networkx, numba, pandas

## Running Tests

```bash
cd mmhc
pytest --cov=src/mmhc tests/ -v
```

## Linting

```bash
cd mmhc
ruff check src/ tests/
```

## Architecture

All source code lives in `mmhc/src/mmhc/`:

- **`__init__.py`** — Public API: `MMHC`, `mmhc()`, `make_student`, `make_rainy`, `MMHCConfig`, `MMHCResult`
- **`_types.py`** — `MMHCConfig` and `MMHCResult` dataclasses
- **`independence.py`** — G-test conditional independence using vectorized contingency tables (`np.ravel_multi_index` + `np.bincount`). No conditioning variable limit.
- **`scoring.py`** — BDeu scoring with Numba-accelerated lgamma loops
- **`mmpc.py`** — MMPC skeleton learning: Forward/Backward phases with `itertools.combinations` for conditioning subsets
- **`hillclimb.py`** — Hill-climbing edge direction with cycle detection via DFS
- **`graph_utils.py`** — DAG validation, cycle detection, adjacency helpers
- **`mmhc.py`** — Top-level orchestrator with scikit-learn-style `fit()` API
- **`datasets.py`** — `make_student()` and `make_rainy()` synthetic data generators
- **`plotting.py`** — NetworkX/matplotlib DAG visualization

### Key Parameters (configurable via MMHCConfig)
- `alpha = 0.05` — significance level for conditional independence tests
- `eta = 1.0` — equivalent sample size for BDeu score
- `max_iterations = 100` — maximum hill-climbing iterations
- `early_stop_rounds = 5` — stop after N non-improving rounds
- `random_seed = None` — seed for reproducibility
