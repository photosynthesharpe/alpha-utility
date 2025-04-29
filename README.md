# alpha-utility
Code for wrangling inputs and outputs of AlphaFold3.

## Formatting
`alpha_format.py` handles wrangling for 3 kinds of inputs: `folding_only`, `one_v_many`, and `many_v_many`.
* `folding_only`: creates one json per protein in a given fasta file.
* `one_v_many`: given two fasta files, creates one json for every pair of proteins across the two files. "One" can actually be an arbitrary number.
*  `many_v_many`: creates one json for every possible pairwise combination of proteins in a single fasta file.

All of these options support the addition of ligands and cofactors; a json will be greated for every combination of (protein/protein pair, ligand, cofactor). Note that this is a combinatoric operation and can therefore get prohibitively large very quickly; implementing some parallelization here would likely be useful, I think it's as straightforward as using Python's `pool` module.

`alpha_format.py` has also been extended to handle multi-subunit proteins. In order to use this functionality, a json specifying each complex and its subunits is required. An example of this json for rubisco:

```
{
    "rubisco": {
        "O03042": 8,
        "P10795": 8
    }
}
```

The outer keys, like "rubisco" can be any arbitrary ID string without spaces. The inner keys must correspond to keys that are present in one of the fasta files provided to `alpha_format.py`. The program assumes that all subunits will be found either in the main fasta file, OR the anchor fasta file, but not both. If not all subunits are present in the dataset, the code will run without raising an error, and will exclude subunits that are not present in the dataset. Additional multimers for the same dataset can be specified as additional key:value pairs.
