#!/bin/bash

dataDir=/home/data/madlab/McMakin_EMUR01/dset
parDir=/scratch/madlab/nate_ppi

# copy data
cd $dataDir
for i in s*; do
    if [ -f ${i}/ses-S2/func/${i}_ses-S2_task-rest_run-1_bold.nii.gz ]; then

        anatDir=${parDir}/dset/${i}/ses-S2/anat
        funcDir=${parDir}/dset/${i}/ses-S2/func
        derivDir=${parDir}/derivatives/${i}/ses-S2
        mkdir -p $anatDir $funcDir $derivDir

        cp ${i}/ses-S2/func/${i}_ses-S2_task-rest_run-*_bold.* $funcDir
        cp ${i}/ses-S1/anat/${i}_ses-S1_run-1_T1w.* $anatDir
    fi
done

# rename
cd ${parDir}/dset
for i in s*; do
    mv ${i}/ses-S2/anat/${i}_ses-S1_run-1_T1w.nii.gz ${i}/ses-S2/anat/${i}_ses-S2_T1w.nii.gz
done
