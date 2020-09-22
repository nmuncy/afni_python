#!/bin/bash

codeDir=~/compute/afni_python
parDir=/scratch/madlab/nate_vCAT
workDir=${parDir}/derivatives
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/TS0_${time}

mkdir -p $outDir

refDir=/home/data/madlab/Mattfeld_vCAT/sourcedata
subjList=(`ls ${refDir}/*tar.gz`)

for i in ${subjList[@]}; do

    file=${i##*/}

    sbatch \
    -o ${outDir}/output_TS0_${file%%.*}.txt \
    -e ${outDir}/error_TS0_${file%%.*}.txt \
    ${codeDir}/step0_wrap_py.sh $file $codeDir

    sleep 1
done

