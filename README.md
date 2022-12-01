# Efficient Enumeration of Markov Equivalent DAGs

This repository contains the implementation of the algorithms presented in the
paper "Efficient Enumeration of Markov Equivalent DAGs" by Marcel Wienöbst,
Malte Luttermann, Max Bannach, and Maciej Liśkiewicz (AAAI 2023).

## Graph Instance Generation
The following generation approaches are supported and have been applied to
generate the graph instances for the experiments of the paper:

- `julia generate.jl cc` to generate undirected chordal graphs.
The instances are stored in the directory `/instances/cc/`.
- `julia generate.jl cpdag` to generate CPDAGs.
The instances are stored in the directory `/instances/cpdag/`.
- `julia generate.jl pdag` to generate PDAGs.
The instances are stored in the directory `/instances/pdag/`.

> Graphs are represented as `SimpleDiGraph`s from the Graphs package, i.e.,
> undirected edges `a - b` are encoded as two directed edges `a -> b` and
> `b -> a`.

## Enumerate all Markov Equivalent DAGs in a Markov Equivalence Class
To enumerate all Markov equivalent DAGs for an input graph, there
are multiple algorithms available:

- `chickering_enumerate` enumerates all DAGs of an input graph based on
the transformational characterization by Chickering (1995).
- `dfs_enumerate` enumerates all DAGs of an input graph such that two
consecutive DAGs have structural Hamming distance less or equal than three.
- `enumerate_meek` enumerates all DAGs of an input graph based on the rules
proposed by Meek (1995).
- `cpdag_enumerate` (for CPDAGs) and `pdag_enumerate` (for PDAGs) enumerate
all DAGs with linear-time delay.

To benchmark all algorithms against each other in every setting, run
`julia run_experiments.jl all`.
To run experiments only for a specific set of graph instances, one
can also execute `julia run_experiments.jl [TYPE]` with `[TYPE]` being
one of `cc`, `cpdag`, `pdag`.