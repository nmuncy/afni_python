#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "TS0"   # job name
#SBATCH --partition centos7_IB_44C_512G
#SBATCH --account iacc_madlab

h_arg=$1
h_dir=$2

module load python-3.7.0-gcc-8.2.0-joh2xyk 
python ${h_dir}/step0_dcm2nii.py $h_arg