#!/bin/bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "TS3"   # job name
#SBATCH --partition centos7_IB_44C_512G
#SBATCH --account iacc_madlab

subj=$1
sess=$2
phase=$3
decon=$4
deriv_dir=$5
code_dir=$6

module load python-3.7.0-gcc-8.2.0-joh2xyk
python ${code_dir}/step3_decon.py $subj $sess $phase $decon $deriv_dir