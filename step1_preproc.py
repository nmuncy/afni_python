# %%
# --- Notes
#
# 1) only written to accept one experiment phase (e.g. Study, Test) per session, assumes
#       all epi scans in session pertain to phase


import json
import os
import sys
import subprocess
import fnmatch
import math


# %%
# --- Step 0: Set up

test_mode = True

if test_mode:
    subj = "sub-4001"
    sess = "ses-S1"
    phase = "Study"
else:
    subj = str(sys.argv[1])
    sess = str(sys.argv[2])
    phase = str(sys.argv[3])

data_dir = os.path.join("/home/data/madlab/McMakin_EMUR01/dset", subj, sess)
work_dir = os.path.join("/scratch/madlab/chen_analysis/derivatives", subj, sess)
atlas_dir = "TODO"

if not os.path.exists(work_dir):
    os.makedirs(work_dir)


# Submit jobs to slurm
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_sub, h_ses):
    full_name = "TP1_" + h_sub + "-" + h_ses
    sbatch_job = "sbatch \
        -J TP1 -t {}:00:00 --mem={}000 --ntasks-per-node={} \
        -p centos7_IB_44C_512G  -o {}.out -e {}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap='module load afni-20.2.06 \n {}'".format(
        wall_hours, mem_gig, num_proc, full_name, full_name, command
    )
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id, error = sbatch_response.communicate()
    return job_id


# Submit afni subprocess
def func_afni(cmd):
    afni_job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    afni_out, afni_err = afni_job.communicate()
    return afni_err


# %%
# --- Step 1: Copy data into work_dir
#
# Get func, anat, fmap data. Rename appropriately.
#
# To account for different num of fmaps, will
#   produce AP, PA fmap per run

# struct
struct_nii = os.path.join(data_dir, "anat", "{}_{}_run-1_T1w.nii.gz".format(subj, sess))
struct_raw = os.path.join(work_dir, "struct+orig")
if not os.path.exists(struct_raw + ".HEAD"):
    h_cmd = "3dcopy {} {}".format(struct_nii, struct_raw)
    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess)
    else:
        func_afni(h_cmd)

# epi - only study, not rest
epi_list = [
    epi
    for epi in os.listdir(os.path.join(data_dir, "func"))
    if fnmatch.fnmatch(epi, "*study*.nii.gz")
]

# start dict for header info
epi_dict = {}
for i in range(len(epi_list)):

    run_num = 1 + i
    str_raw = "run-{}_{}+orig".format(run_num, phase)
    epi_raw = os.path.join(work_dir, str_raw)
    epi_nii = os.path.join(data_dir, "func", epi_list[i])

    with open(epi_nii.split(".")[0] + ".json") as j:
        epi_dict[str_raw.split("+")[0]] = json.load(j)

    if not os.path.exists(epi_raw + ".HEAD"):
        h_cmd = "3dcopy {} {}".format(epi_nii, epi_raw)
        if test_mode:
            func_sbatch(h_cmd, 1, 1, 1, subj, sess)
        else:
            func_afni(h_cmd)

# fmap
json_list = [
    x
    for x in os.listdir(os.path.join(data_dir, "fmap"))
    if fnmatch.fnmatch(x, "*.json")
]

fmap_list = []
for i in json_list:
    with open(os.path.join(data_dir, "fmap", i)) as j:
        h_json = json.load(j)
        for k in epi_list:
            h_epi = os.path.join(sess, "func", k)
            if h_epi in h_json["IntendedFor"]:
                fmap_list.append(i.split(".")[0] + ".nii.gz")

# copy fmap function
def func_fmap(h_file, h_run):

    fmap_nii = os.path.join(data_dir, "fmap", h_file)
    h_dir = h_file.split("-")[4].lstrip().split("_")[0]
    enc_dir = "Forward" if h_dir == "AP" else "Reverse"
    fmap_raw = os.path.join(work_dir, "blip_run-{}_{}+orig".format(h_run, enc_dir))

    if not os.path.exists(fmap_raw + ".HEAD"):
        h_cmd = "3dcopy {} {}".format(fmap_nii, fmap_raw)
        if test_mode:
            func_sbatch(h_cmd, 1, 1, 1, subj, sess)
        else:
            func_afni(h_cmd)


# TODO: first half of conditional is untested

if len(fmap_list) == 2:
    for i in range(1, 3):
        for j in fmap_list:
            func_fmap(j, i)
else:
    h_half = len(fmap_list) // 2
    fmap_list_A = sorted(fmap_list[:h_half])
    fmap_list_B = sorted(fmap_list[h_half:])

    count = 1
    for i, j in zip(fmap_list_A, fmap_list_B):
        func_fmap(i, count)
        func_fmap(j, count)
        count += 1


# %%
# --- Step 2: Detect outliers voxels, blip correct
#
# num_tr is hardcoded (305), this doesn't exist in json file

for i in epi_dict.keys():

    pol = 1 + math.ceil((epi_dict[i]["RepetitionTime"] * 305) / 150)
    h_fileA = os.path.join(work_dir, "outcount." + i + ".1D")
    h_fileB = os.path.join(work_dir, "out.cen." + i + ".1D")

    h_cmd = """
        3dToutcount -automask -fraction -polort {} -legendre {}+orig > {}
        1deval -a {} -expr "1-step(a-0.1)" > {}
    """.format(
        pol,
        os.path.join(work_dir, i),
        h_fileA,
        h_fileA,
        h_fileB,
    )

    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess)
    else:
        func_afni(h_cmd)


# ### Update - do blip here
# # Ripped from afni_proc.py

# if [ ! -f ${block[0]}_blip+orig.HEAD ]; then

# 	# create median datasets from forward and reverse time series
# 	3dTstat -median -prefix rm.blip.med.fwd blip_Forward+orig
# 	3dTstat -median -prefix rm.blip.med.rev blip_Reverse+orig

# 	# automask the median datasets
# 	3dAutomask -apply_prefix rm.blip.med.masked.fwd rm.blip.med.fwd+orig
# 	3dAutomask -apply_prefix rm.blip.med.masked.rev rm.blip.med.rev+orig

# 	# compute the midpoint warp between the median datasets
# 	3dQwarp -plusminus -pmNAMES Rev For                           \
# 	    -pblur 0.05 0.05 -blur -1 -1                          \
# 	    -noweight -minpatch 9                                 \
# 	    -source rm.blip.med.masked.rev+orig                   \
# 	    -base   rm.blip.med.masked.fwd+orig                   \
# 	    -prefix blip_warp

# 	# # warp median datasets (forward and each masked) for QC checks
# 	# # (and preserve obliquity)
# 	# 3dNwarpApply -quintic -nwarp blip_warp_For_WARP+orig          \
# 	# 	-source rm.blip.med.fwd+orig                     \
# 	# 	-prefix blip_med_for

# 	# 3drefit -atrcopy blip_forward+orig IJK_TO_DICOM_REAL          \
# 	# 	blip_med_for+orig

# 	# 3dNwarpApply -quintic -nwarp blip_warp_For_WARP+orig          \
# 	# 	-source rm.blip.med.masked.fwd+orig              \
# 	# 	-prefix blip_med_for_masked

# 	# 3drefit -atrcopy blip_forward+orig IJK_TO_DICOM_REAL          \
# 	# 	blip_med_for_masked+orig

# 	# 3dNwarpApply -quintic -nwarp blip_warp_Rev_WARP+orig          \
# 	# 	-source rm.blip.med.masked.rev+orig              \
# 	# 	-prefix blip_med_rev_masked

# 	# 3drefit -atrcopy blip_reverse+orig IJK_TO_DICOM_REAL          \
# 	# 	blip_med_rev_masked+orig

# 	# warp EPI time series data
# 	for i in ${block[@]}; do
# 	    3dNwarpApply -quintic -nwarp blip_warp_For_WARP+orig      \
# 			-source ${i}+orig           \
# 			-prefix ${i}_blip

# 	    3drefit -atrcopy blip_forward+orig IJK_TO_DICOM_REAL      \
# 			${i}_blip+orig
# 	done
# fi

# %%
