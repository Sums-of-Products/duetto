# Duetto
Duetto implements **Markov chain Monte Carlo (MCMC)** moves based on **duet sampling** for sampling **directed acyclic graphs (DAGs)** from their posterior distribution. 

## Algorithm Overview
Duetto implements three move types

- **Unconstrained Sequential Resampling (USR):** Sequentially resamples parent sets following the order given by a node permutation.  
- **Edge-constrained move (EC):** Adds, removes, or reverses an edge while resampling parent sets at its endpoints.
- **Ancestor-constrained move (AC):** Adds, removes, or reverses an ancestor relationship while resampling parent sets at its endpoints.
- **[Blocked Gibbs sampling](https://www.jmlr.org/papers/v17/14-486.html) (BG)** (no duet sampling here) 

A **score-pruning technique** can optionally be enabled to discard low-probability parent sets during local score computation, improving scalability and memory efficiency.

## Compilation and Execution

Gibby is implemented in **C++** and can be compiled with **g++** as follows:
