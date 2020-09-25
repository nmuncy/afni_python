#!/bin/bash

codeDir=~/compute/afni_python
workDir=/scratch/madlab/nate_vCAT/derivatives
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/TS3_${time}

mkdir -p $outDir

sessList=(ses-S1)
phaseList=(vCAT)
deconType=2GAM

cd $workDir
for i in s*; do
    for j in ${sessList[@]}; do
        for k in ${phaseList[@]}; do

            sbatch \
            -o ${outDir}/output_TS3_${i}_${j}_${k}.txt \
            -e ${outDir}/error_TS3_${i}_${j}_${k}.txt \
            ${codeDir}/step1_wrap_py.sh $i $j $k $deconType $workDir $codeDir

            sleep 1
        done
    done
done

