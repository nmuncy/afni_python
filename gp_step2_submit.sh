#!/bin/bash

module load r-3.5.1-gcc-8.2.0-djzshna

behavDir=/home/data/madlab/Mattfeld_vCAT/behav
derivDir=/scratch/madlab/nate_vCAT/derivatives
codeDir=~/compute/afni_python

# # For testing
# behavDir=/Users/nmuncy/Projects/learn_mvpa/vCAT_data
# derivDir=/Users/nmuncy/Projects/afni_python
# codeDir=$derivDir

numRuns=2

cd $behavDir
for i in vCAT*; do

	subj=sub-${i#*_}
	dataDir=${behavDir}/$i
	outDir=${derivDir}/${subj}/ses-S1

	if [ -d $outDir ]; then
		Rscript ${codeDir}/gp_step2_timingFiles_localizer.R $dataDir $outDir $i $numRuns
	fi
done
