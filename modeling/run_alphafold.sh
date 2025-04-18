#!/bin/bash

program_name=$0

function Help() {
    echo
    echo "Run AlphaFold3 over a directory of jsons using MSU's HPCC. "
    echo
    echo "Syntax: $program_name [-h] [-p parameter_dir] [-j json_dir] [-o outdir]"
    echo
    echo "Options:"
    echo "  -p parameter_dir specify path to your copy of the AlphaFold3 "
    echo "                       parameters. This should be the directory "
    echo "                       that contains the af3.bin.zst file from "
    echo "                       Google. Provide an absolute path."
    echo "  -j json_dir      specify the directory containing jsons to run"
    echo "  -o outdir        specify path to save output"
}

# Define accepted args
while getopts ":hp:j:o:" option; do
    case $option in
      h)
         Help
         exit;;
      p)
        parameter_dir=$OPTARG;;
      j) 
        json_dir=$OPTARG;;
      o)
        outdir=$OPTARG;;
      :)
        echo "Option -${option} requires an argument."
        exit 1;;
      ?) # Invalid option
        echo "Error: Invalid option -${option}"
        exit 1;;
    esac
done

# First, purge and load modules
module purge
module load AlphaFold3/3.0.1-foss-2023a-CUDA-12.4.0 ## Could change

# Then define the parameter and database global variables
export MODEL_WEIGHTS=$parameter_dir
export DATABASES="/mnt/research/common-data/alphafold/database_3/"

# Run for all jsons in the directory
command="run_alphafold.py \
  --input_dir=$json_dir \
  --model_dir=${MODEL_WEIGHTS}
  --db_dir=${DATABASES}
  --output_dir=$outdir"
echo $command