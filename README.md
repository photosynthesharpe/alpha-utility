# alpha-utility
Code for wrangling inputs and outputs of AlphaFold3.

There are 3 major components of this pipeline: **formatting**, **modeling**, and **analysis**.

## Formatting

AlphaFold3 requires a specific format of `.json` file as an input for each fold. More details on the format can be found [here](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md); however, the script `alpha_format.py` in the `formatting` directory is designed to automatically generate input json files for a variety of use cases.

`alpha_format.py` handles wrangling for 3 kinds of inputs: `folding_only`, `one_v_many`, and `many_v_many`.
* `folding_only`: creates one json per protein in a given fasta file.
* `one_v_many`: given two fasta files, creates one json for every pair of proteins across the two files. "One" can actually be an arbitrary number.
*  `many_v_many`: creates one json for every possible pairwise combination of proteins in a single fasta file.

All of these options support the addition of ligands and cofactors; a json will be greated for every combination of (protein/protein pair, ligand, cofactor). Note that this is a combinatoric operation and can therefore get prohibitively large very quickly; implementing some parallelization here would likely be useful, I think it's as straightforward as using Python's `pool` module.

The code has recently (11/21/2025) been updated to support the addition of multiple ligands and cofactors per fold. If you'd like all provided ligands and/or cofactors to be included in each json, use the parameters `--group_ligands` and ``group_cofactors`.

A basic example of how to run `alpha_format.py` for the `one_v_many` case:

```
python alpha_format.py /path/to/proteins.fasta /path/to/save/output/ <file_prefix> -anchor_fasta /path/to/anchor.fasta -protein_comparison_type one_v_many
```

`alpha_format.py` has also been extended to handle multi-subunit proteins. In order to use this functionality, a json specifying each complex and its subunits is required. An example of this json for rubisco:

```
{
    "rubisco": {
        "O03042": 8,
        "P10795": 8
    }
}
```

The outer keys, like "rubisco" can be any arbitrary ID string without spaces. The inner keys must correspond to keys that are present in one of the fasta files provided to `alpha_format.py`. The keys are the string that appears after the `>` in the fasta file headers, but before the first space. The program assumes that all subunits will be found either in the main fasta file, OR the anchor fasta file, but not both. If not all subunits are present in the dataset, the code will run without raising an error, and will exclude subunits that are not present in the dataset. Additional multimers for the same dataset can be specified as additional key:value pairs. To facilitate making a json, I'd reccomend typing or pasting your multimer string into a website like [JSONLint](https://jsonlint.com/) to catch any formatting errors.

An example of how to run the `one_v_many` case with multimers:
```
python alpha_format.py /path/to/proteins.fasta /path/to/save/output/ 16June2025_random -anchor_fasta /path/to/anchor.fasta -protein_comparison_type one_v_many -multi_subunit_proteins /path/to/multimers.json
```

Note that if you'd like to fold more than two proteins against one another, like several monomers plus a complex, you can create artificial "complexes" in the input complex json. This will have the same final effect as suplying multiple unique monomers to the same json. The complex input file will have a key that you make up to name that artificial complex, but otherwise will look functionally identical to a real complex like the rubisco example above.

## Modeling

The modeling directory contains templates for job submission scripts for Michigan State's High-Percormance Computing Center, though they should also work on other slurm systems. In order to shorten the total time our jobs wait in the queue, we will separate the multiple sequence alignment (MSA) part of AlphaFold (the data pipeline) and the prediction part (inference). MSA is the time-consuming part, but does not use GPU, so we can make a job submission that gets out of the queue much faster by doing the time-intensive portion without requesting valuable resources like GPU. We then run the inference component later, so we can be in the short queue for GPU and spend less time waiting for resources.

Scripts in this directory are: 
* `00_run_sequence_alignment_EXAMPLE.sb`: Do the CPU time-intensive portion
* `01_run_inference_EXAMPLE.sb`: Do the GPU inference portion

### Customizing templates
To use these scripts, copy to another (non-tracked) directory, and replace anything within `<>` with the relevant path or information. Specifically;

In `00_run_sequence_alignment_EXAMPLE.sb`, the items to replace for a Walker Lab user are:
* `<NUM_FILES>` -- see below for how to obtain the number that goes here
* `<DATA_DIR_IN_QUOTES>` -- replace this with the file path that contains the outputs to the formatting script (the json files produced by `alpha_format.py`
Other MSU users will also need to change the path to the `MODEL_WEIGHTS` parameter after downloading your own copy of the AlphaFold model weights [here](MODEL_WEIGHTS).

For non-MSU users on slurm systems, there are some additional fields you may need to change:
* `DATABASES` -- MSU's HPC has a global installation of the large databases required to run the MSA step; on other systems you may have to create your own installation or provide a path to a different global installation.
* `module load AlphaFold3/3.0.1-foss-2023a-CUDA-12.4.0` -- you'll need to change this to whatever the version of AlphaFold3 is available on your system.

In `01_run_inference_EXAMPLE.sb`, you'll need to replace all the same parameters; however, see note below about `<NUM_FILES>`.

### Array job-specific parameters (`<NUM_FILES>`)
Because these are array jobs, you'll need to know how many files are in your input directory; unfortunately this can't be done automatically in the job script, because slurm stops reading `#SBATCH` directives once any script command has been executed. So you'll need to run the following on the command line in the directory you used as `<DATA_DIR_IN_QUOTES>` before submitting any jobs (`cd` to `<DATA_DIR_IN_QUOTES>` before running):
```
ls ./*.json > file_list
wc -l file_list
```
Then put that number in place of `<NUM_FILES>` in the job script. The job script also reads the `file_list` file, so make sure to actually produce it, rather than piping `ls` to `wc`.

If the number of files is greater than 1000, you'll run into issues on job submission, as 1000 is the upper limit of jobs a single user can have in the queue. To deal with this, use the script `split_folds.py`. For the above example, the command would be: `python split_folds.py /<path_to_my_dir>/file_list`. This will create a set of `file_list` files with 999 jsons in each; you would then submit the `00_run_sequence_alignment_EXAMPLE.sb` (and then the inference one after) for each of the generated subset files.

The process of determining how many array jobs to submit and what they are needs to be repeated before running `01_run_inference_EXAMPLE.sb`. Because the input json names aren't necessarily the same as the output directory/file names, we regenerate the input list for each of the directories with the following from within the data directory:

```
ls -d */* > data_file_list
wc -l data_file_list
```
As before, if there are more than 1000 files, you'll need to split them with the `split_folds.py` script.

Note that the inference outputs will appear in new directories. These will have the same base name as the ones that the data pipeline output, but they have a timestamp added to them in the format `_YYYYMMDD_<time>`. This is [intentional](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md#output-directory-structure) to avoid overwriting data; just know that the final inference results you want are the ones from the timestamped directories.

## Analysis
The idea behind the scripts in the `analysis` directory is to use metrics generated by AlphaFold, plus RMSD, to identify high-likelihood protein interactions. For analysis including rubisco, the script `analysis/symmetricalRMSD.py` handles the D4 symmetry to prevent artificial inflation of RMSD. Approaches in the `analysis` directory should be used with caution -- please see the [CoIP notebooks](https://github.com/photosynthesharpe/alpha-utility/blob/main/notebooks/coIP_vs_predicted.ipynb) for an analysis of these metrics.
