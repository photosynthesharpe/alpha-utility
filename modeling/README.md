# Job submission scripts
The scripts in this directory are templates for how to write job submission scripts for Michigan State's High-Percormance Computing Center. To use them, copy to another (non-tracked) directory, and replace anything within `<>` with the relevant path or information.

In order to shorten the total time our jobs wait in the queue, we will separate the multiple sequence alignment (MSA) part of AlphaFold (the data pipeline) and the prediction part (inference). MSA is the time-consuming part, but does not use GPU, so we can make a job submission that gets out of the queue much faster by doing the time-intensive portion without requesting valuable resources like GPU. We then run the inference component later, so we can be in the short queue for GPU and spend less time waiting for resources.

Scripts in this directory are: 
* `00_run_sequence_alignment_EXAMPLE.sb`: Do the CPU time-intensive portion
* `01_run_inference_EXAMPLE.sb`: Do the GPU inference portion

Because these are array jobs, you'll need to know how many files are in your input directory; unfortunately this can't be done automatically in the job script, because slurm stops reading `#SBATCH` directives once any script command has been executed. So you'll need to run the following on the command line in your directory of interest before submitting any jobs:
```
ls ./*.json > file_list
wc -l file_list
```
Then put that number in place of `<NUM_FILES>` in the job script. The job script also reads the `file_list` file, so make sure to actually produce it, rather than piping `ls` to `wc`.

If the number of files is greater than 1000, you'll run into issues on job submission, as 1000 is the upper limit of jobs a single user can have in the queue. To deal with this, use the script `split_folds.py`. This will create a set of `file_list` files with 999 jsons in each; you would then submit the `00_run_sequence_alignment_EXAMPLE.sb` (and then the inference one after) for each of the generated subset files.
