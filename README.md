# Duetto
Duetto implements **Markov chain Monte Carlo (MCMC)** moves based on **duet sampling** to draw from the posterior distribution of **directed acyclic graphs (DAGs)**. The underlying framework and moves are introduced in **Generalizing Metropolis—Hastings for Sampling DAGs**.

## Algorithm Overview
Duetto implements four moves:

- **Unconstrained Sequential Resampling (USR):** Sequentially resamples parent sets following the order given by a node permutation.  
- **Edge-constrained move (EC):** Adds, removes, or reverses an edge while resampling parent sets at its endpoints.
- **Ancestor-constrained move (AC):** Adds, removes, or reverses an ancestor relationship while resampling parent sets at its endpoints.
- **[Blocked Gibbs sampling](https://www.jmlr.org/papers/v17/14-486.html) (BG):** standard Gibbs sampling (no duet sampling)

The **score-pruning** technique, introduced in **[Scaling Up Bayesian DAG Sampling](https://proceedings.mlr.press/v337/nikzad26a.html)** (Nikzad et al.) can be enabled to discard low-probability parent sets in the preprocessing phase, improving scalability and memory efficiency.

## Compilation and Execution

Gibby is implemented in **C++** and can be compiled with **g++** as follows:
