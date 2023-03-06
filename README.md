# Heuristics for Cycle Packing of Adjacency Graphs for Genomes with Repeated Genes

This repository has the implementation of two heuristics for the Adjacency Graph Packing problem described in the paper ``Heuristics for Cycle Packing of Adjacency Graphs for Genomes with Repeated Genes". One is called Random Packings and tries to solve the problem by randomly generating cycle packings of the Adjacency Graph. The other is based on the Genetic Algorithm metaheuristic.

The repository also presents the implementation of multiple algorithms to approximate rearrangement distances between strings. The following algorithms are implemented:
- A 3 asymptotic approximation algorithm for Sorting by Intergenic Transposition.
- An 4-approximation algorithm for Sorting by Intergenic Reversal [[1]](#1).
- An 4.5-approximation algorithm for Sorting by Intergenic Reversal and Transposition [[1]](#1).
- An 2-approximation for Reversal and Indel Distance in signed strings. This algorithm uses an implementation of an exact algorithm for Sorting by Reversal [[2]](#2).
- An heuristic for Reversal, Transposition and Indels Distance in signed strings. This algorithm uses an implementation of a 2-approximation algorithm for Sorting by Reversal and Transposition [[3]](#3).

The algorithms can be applied to genomes with multiple genes, in that case random mappings of the genomes into permutations will be used to produce the distances.

## Usage

Compile the code by running `make` and see the running options with `./dec --help` (for the cycle packing) or `./dist --help` (for the rearrangement distances). To use other algorithm when calculating the distances include a executable in the `external` folder and pass its name as the algorithm parameter. In that case, the executable should receive one instance as command line arguments (four coma separated lists) and produce the distance in the standard output.

## References

<a id="1">[1]</a> 
Klairton Lima Brito, Géraldine Jean, Guillaume Fertin, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. Sorting by genome rearrangements on both gene order and intergenic sizes. Journal of Computational Biology 2020;27(2):156–74.

<a id="2">[2]</a> 
Anne Bergeron, Julia Mixtacki and Jens Stoye. “Reversal Distance without Hurdles and Fortresses”. Combinatorial Pattern Matching. Berlin, Heidelberg: Springer.

<a id="3">[3]</a> 
Maria Emília M.T. Walter, Zanoni Dias and João Meidanis. “Reversal and transposition distance of linear chromosomes”. Proceedings of the String Processing and Information Retrieval: A South American Symposium (SPIRE’1998). Los Alamitos, CA, USA, 1998, pp. 96–102.
