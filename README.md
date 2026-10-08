# Duetto

Duetto implements **Markov chain Monte Carlo (MCMC)** moves based on **duet sampling** to draw from the posterior distribution of **directed acyclic graphs (DAGs)**. The underlying framework and moves are introduced in **Generalizing Metropolis—Hastings for Sampling DAGs**.

## Algorithm Overview

Duetto implements four moves:

- **Unconstrained Sequential Resampling (USR):** Sequentially resamples parent sets following the order given by a node permutation.  
- **Edge-constrained move (EC):** Adds, removes, or reverses an edge while resampling parent sets at its endpoints.
- **Ancestor-constrained move (AC):** Adds, removes, or reverses an ancestor relationship while resampling parent sets at its endpoints.
- **[Blocked Gibbs sampling](https://www.jmlr.org/papers/v17/14-486.html) (BG):** jointly samples a block of parent sets from the full conditional DAG posterior (no duet sampling).

The **score-pruning** technique, introduced in **[Scaling Up Bayesian DAG Sampling](https://proceedings.mlr.press/v337/nikzad26a.html)** (Nikzad et al., 2026) can be enabled to discard low-probability parent sets in the preprocessing phase, improving scalability and memory efficiency.

## Compilation and Execution

Duetto is implemented in **C++** and can be compiled with **g++** as follows:

```bash
g++ -std=c++17 -march=native -O3 -o duetto duetto.cpp
```

The only mandatory parameter is the data file, which should contain data points as rows and features as columns. Alternatively, the data can be provided with the opposite orientation by using the `.csv` extension. The `.csv` extension is required for the file to be correctly interpreted with this orientation.

Another option is to provide a local scores file in `.jkl` format. The `.jkl` extension is required for the file to be correctly detected as a score file.

For example, run:

```bash
./gibby  ./data/asia10k.dat 
```

The datasets and score files used in the experiments are located in the `data` folder.

Duetto computes scores only from discrete data using BDeu (ESS=1). We recommend using [GOBNILP](https://benchpressdocs.readthedocs.io/en/latest/structure_learning_algorithms/gobnilp.html) to generate local score files for continuous data. 

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-K` | Top candidate parents per node | `0` (unrestricted) |
| `-m` | Maximum in-degree per node | `4 4` The first value is for the parent sets containing the top K nodes, the second is for the others |
| `-p` | Structure prior: `0` uniform, `1` fair, `2` fair+, `-w` edge(w) | `1` |
| `-P` | Score pruning mode: `0` none, `1` top-down, `2` bottom-up; if preceded by digit `k>0`, only up to `min{k, maxind}` parents | `0` |
| `-a` | Number of significant bits in approximations | `15` |
| `-R` | Seed for random number generator | random |
| `-i` | Number of main iterations | `10 000` |
| `-burnin` | Number of burn-in iterations | `i / 10` |
| `-q` | Block size | $\min(\lfloor\sqrt{\text{n.nodes}}\rfloor, 10)$ |
| `-USRrev` | Number of USRrev moves per iteration | `100` |
| `-USRrnd` | Number of USRrnd moves per iteration | `50` |
| `-USRid` | Number of USRrev moves per iteration | `0` |
| `-EC` | Number of EC moves per iteration | `25` |
| `-AC` | Number of AC moves per iteration | `25` |
| `-BG` | Number of BG moves per iteration | `0` |
| `-iPPE` | Generates a file showing edge posterior probabilities every iPPE `<n>` iterations | none |


## Output Files

After execution, duetto produces the following outputs:

- **Sampled DAGs** — contains the sampled DAGs. Each DAG is listed as node–parent relationships in the format `node: parent1,parent2,...`. Nodes with empty parent sets are omitted.
- **Sampled DAG scores** — contains the score of the sampled DAGs.
- **Edge probability matrix** — contains the posterior probabilities of the edges, where rows represent parents and columns represent children.
- **(Optional) Edge probability matrix per iPPE** — if -iPPE `<n>` is specified, a file showing edge posterior probabilities every n iterations is generated.


