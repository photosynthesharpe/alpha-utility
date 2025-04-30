# Job submission scripts
The scripts in this directory are templates for how to write job submission scripts for Michigan State's High-Percormance Computing Center. To use them, copy to another (non-tracked) directory, and replace anything within `<>` with the relevant path or information.

Scripts in this directory are: 
* `EXAMPLE_run_alphafold.sb`: uses one GPU to run AlphaFold3 over all json files in an input directory without parallelization.
* `EXAMPLE_run_alphafold_PARALLEL.sb`: uses a job array to run all json files in an input directory in parallel.

To run the array job, you'll need to know how many files are in your input directory; unfortunately this can't be done automatically in the job script, because slurm stops reading `#SBATCH` directives once any script command has been executed. So you'll need to run the following on the command line in your directory of interest:
```
ls ./*.json | wc -l
```
Then put that number in place of `<NUM_FILES>` in the job script.