# alpha-utility
Code for wrangling inputs and outputs of AlphaFold3.

## Formatting
`alpha_format.py` handles wrangling for 3 kinds of inputs: `folding_only`, `one_v_many`, and `many_v_many`.
* `folding_only`: creates one json per protein in a given fasta file.
* `one_v_many`: given two fasta files, creates one json for every pair of proteins across the two files. "One" can actually be an arbitrary number.
*  `many_v_many`: creates one json for every possible pairwise combination of proteins in a single fasta file.
All of these options support the addition of ligands and cofactors; a json will be greated for every combination of (protein/protein pair, ligand, cofactor). Note that this is a combinatoric operation and can therefore get prohibitively large very quickly; implementing some parallelization here would likely be useful, I think it's as straightforward as using Python's `pool` module.
