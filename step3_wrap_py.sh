#!/bin/bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "TS3"   # job name
#SBATCH --partition centos7_IB_44C_512G
#SBATCH --account iacc_madlab

h_subj=$1
h_sess=$2
h_phase=$3
# h_dir=$4
h_dir=~/compute/afni_python

module load python-3.7.0-gcc-8.2.0-joh2xyk
python ${h_dir}/step3_decon.py