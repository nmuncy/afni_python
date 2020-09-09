# %%

import json
import os
import sys
import subprocess
import fnmatch


# %%
# Receive/Set orienting vars

# subj = str(sys.argv[1])
# sess = str(sys.argv[2])
# phase = str(sys.argv[3])

subj = "sub-4001"
sess = "ses-S1"
phase = "Study"

test_mode = True

data_dir = os.path.join("/home/data/madlab/McMakin_EMUR01/dset", subj, sess)
work_dir = os.path.join("/scratch/madlab/chen_analysis/derivatives", subj, sess)
atlas_dir = "TODO"

if not os.path.exists(work_dir):
    os.makedirs(work_dir)


# %%
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


def func_afni(cmd):
    afni_job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    afni_out, afni_err = afni_job.communicate()
    return afni_err


# %%
# get, rename struct
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

for i in range(len(epi_list)):
    run_num = 1 + i
    epi_raw = os.path.join(work_dir, "run-{}_{}+orig".format(run_num, phase))
    epi_nii = os.path.join(data_dir, "func", epi_list[i])
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

# for simplicity, make AP/PA for each run
if len(fmap_list) == 2:
    for i in range(1, 3):
        for j in fmap_list:

            fmap_nii = os.path.join(data_dir, "fmap", j)
            h_dir = j.split("-")[4].lstrip().split("_")[0]
            enc_dir = "Forward" if h_dir == "AP" else "Reverse"
            fmap_raw = os.path.join(work_dir, "blip_run-{}_{}+orig".format(i, enc_dir))

            if not os.path.exists(fmap_raw + ".HEAD"):
                h_cmd = "3dcopy {} {}".format(fmap_nii, fmap_raw)
                if test_mode:
                    func_sbatch(h_cmd, 1, 1, 1, subj, sess)
                else:
                    func_afni(h_cmd)
else:
    h_half = len(fmap_list) // 2
    fmap_list_A = fmap_list[: len(fmap_list) // 2]  ### Need to sort
    fmap_list_B = fmap_list[len(fmap_list) // 2 :]
    count = 1
    for i, j in fmap_list_A, fmap_list_B:
        print(i, j)
        fmap_nii_i = os.path.join(data_dir, "fmap", i)
        h_dir_i = i.split("-")[4].lstrip().split("_")[0]

# %%
