# MMHC — Max-Min Hill-Climbing Algorithm

A Python implementation of the Max-Min Hill-Climbing (MMHC) Bayesian network structure learning algorithm.

Based on: Tsamardinos, Brown & Aliferis, "The max-min hill-climbing Bayesian network structure learning algorithm", Machine Learning, 2006. DOI: 10.1007/s10994-006-6889-7.

## Installation

```bash
pip install -e ".[dev]"
```

## Quick Start

```python
from mmhc import MMHC, make_student, MMHCConfig

# Generate synthetic data from a known Bayesian network
data = make_student(5000, random_state=42)

# Learn the network structure
model = MMHC(config=MMHCConfig(random_seed=42))
result = model.fit(data)

# Inspect results
print(result.adjacency_matrix)
print(f"BDeu score: {result.score}")
print(f"Converged: {result.converged} in {result.n_iterations} iterations")
```

Or use the convenience function:

```python
from mmhc import mmhc, make_student

data = make_student(5000, random_state=42)
result = mmhc(data)
```

## Algorithm

The MMHC algorithm reconstructs Bayesian Networks from observational data in two phases:

### Phase 1: MMPC (Max-Min Parents and Children)

Builds the undirected skeleton using conditional independence tests (G-test / chi-squared):

$$G = 2 \sum_{i,j} O_{ij} \ln\frac{O_{ij}}{E_{ij}}$$

The forward phase iteratively adds the most dependent variable to the candidate parent-children set. The backward phase removes variables that become conditionally independent given subsets of the current set. A symmetry constraint ensures consistency.

### Phase 2: Hill-Climbing with BDeu Scoring

Directs edges using greedy hill-climbing with BDeu (Bayesian Dirichlet likelihood-equivalence uniform) scoring:

$$\text{BDeu}(X_i, \Pi_i) = \sum_{j=1}^{q_i} \left[ \ln\frac{\Gamma(\eta/q_i)}{\Gamma(N_{ij} + \eta/q_i)} + \sum_{k=1}^{r_i} \ln\frac{\Gamma(N_{ijk} + \eta/(q_i r_i))}{\Gamma(\eta/(q_i r_i))} \right]$$

Edge operations (add, reverse, delete) are applied greedily, with cycle detection to maintain DAG validity. Early stopping triggers after 5 non-improving rounds.

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 0.05 | Significance level for conditional independence tests |
| `eta` | 1.0 | Equivalent sample size for BDeu scoring |
| `max_iterations` | 100 | Maximum hill-climbing iterations |
| `early_stop_rounds` | 5 | Stop after N non-improving rounds |
| `random_seed` | None | Seed for reproducibility |

## API Reference

### Classes

- **`MMHC(config=None)`** — Main class with `fit(data)`, `fit_predict(data)` methods
- **`MMHCConfig(...)`** — Configuration dataclass
- **`MMHCResult`** — Result dataclass with `adjacency_matrix`, `parent_children`, `score`, `node_scores`, `n_iterations`, `converged`, `column_names`

### Functions

- **`mmhc(data, config=None)`** — Convenience function for one-shot usage
- **`make_student(n_samples, random_state)`** — Generate student network data
- **`make_rainy(n_samples, random_state)`** — Generate sprinkler/rain network data
- **`plot_dag(adjacency, labels, ax, title)`** — Visualize the learned DAG

## Use Cases

1. **Medical diagnosis networks**: Learn dependencies between symptoms, diseases, and test results from patient records
2. **Gene regulatory network inference**: Discover gene interaction networks from expression data
3. **Supply chain dependency analysis**: Identify causal relationships between supply chain variables

## Running Tests

```bash
pytest --cov=src/mmhc tests/ -v
```

## Key Improvements Over Original R/C++ Implementation

- **No conditioning variable limit**: Original supported max 3 conditioning variables; Python version handles arbitrary numbers via vectorized contingency tables
- **Cycle detection**: DAG constraint enforced during hill-climbing (original did not check)
- **Reproducibility**: Deterministic with seed (original used `srand(time(NULL))`)
- **Configurable parameters**: All hardcoded values now configurable via `MMHCConfig`
