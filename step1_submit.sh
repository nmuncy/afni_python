#!/bin/bash

codeDir=~/compute/afni_python

parDir=/scratch/madlab/nate_vCAT
dsetDir=${parDir}/dset
workDir=${parDir}/derivatives

slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/TS1_${time}

mkdir -p $outDir

sessList=(ses-S1)
phaseList=(vCAT)

cd $dsetDir
# for i in s*; do
for i in sub-006; do
    for j in ${sessList[@]}; do
        for k in ${phaseList[@]}; do
            if [ ! -f ${workDir}/${i}/${j}/run-1_${k}_scale+tlrc.HEAD ]; then

                sbatch \
                -o ${outDir}/output_TS1_${i}_${j}_${k}.txt \
                -e ${outDir}/error_TS1_${i}_${j}_${k}.txt \
                ${codeDir}/step1_wrap_py.sh $i $j $k $codeDir

                sleep 1
            fi
        done
    done
done

