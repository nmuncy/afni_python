# %%
# --- Notes
#
# 1) only written to accept one experiment phase (e.g. Study, Test) per session, assumes
#       all epi scans in session pertain to phase
#
# TODO:
#   1) Update wait function to wait for job name rather than job number
#   2) Add atlas
#   3) Test various functions/syntax


import json
import os
import sys
import subprocess
import fnmatch
import math
import time


# %%
# --- Step 0: Set up
#
# Set/receive arguments, set variables for data, working, and atlas directories,
#   make working dir if needed, and set up general functions.

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
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_sub, h_ses, h_str):

    full_name = "TP1_" + h_sub + "_" + h_ses + "_" + h_str
    sbatch_job = f"sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap='module load afni-20.2.06 \n {command}'"

    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id, error = sbatch_response.communicate()
    # return job_id

    while_count = 0
    status = False
    while not status:

        check_cmd = "squeue -u $(whoami)"
        sq_check = subprocess.Popen(check_cmd, shell=True, stdout=subprocess.PIPE)
        out_lines, err_lines = sq_check.communicate()
        b_decode = out_lines.decode("utf-8")
        num_lines = len(b_decode.split("\n"))

        if num_lines < 3:
            status = True

        # if not h_str in b_decode:
        #     status = True

        if not status:
            while_count += 1
            print(f"Wait count for sbatch job {h_str}: ", while_count)
            time.sleep(2)
    print("Sbatch jobs finished")


# Submit afni subprocess
# TODO: This function is untested
def func_afni(cmd):
    afni_job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).wait()
    afni_out, afni_err = afni_job.communicate()
    return afni_err


# %%
# --- Step 1: Copy data into work_dir
#
# Get func, anat, fmap data. Rename appropriately.
#
# To account for different num of fmaps, will
#   produce AP, PA fmap per run.

# struct
struct_nii = os.path.join(data_dir, "anat", "{}_{}_run-1_T1w.nii.gz".format(subj, sess))
struct_raw = os.path.join(work_dir, "struct+orig")
if not os.path.exists(struct_raw + ".HEAD"):
    h_cmd = "3dcopy {} {}".format(struct_nii, struct_raw)
    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "struct")
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
            func_sbatch(h_cmd, 1, 1, 1, subj, sess, "epi")
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
            func_sbatch(h_cmd, 1, 1, 1, subj, sess, "fmap")
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
# 1) Determine the proportion of voxels per volume that have outlier signal.
#       Censor volumes that exceed limit.
#
# 2) Correct for signal fallout using fmap. This approach is taken from
#       afni_proc. It uses the fmap to "unwarp" the run epi.
#
# num_tr is hardcoded (305), this doesn't exist in json file

for i in epi_dict.keys():

    # determine polort arg, file strings, find outliers
    pol = 1 + math.ceil((epi_dict[i]["RepetitionTime"] * 305) / 150)
    h_fileA = os.path.join(work_dir, "outcount." + i + ".1D")
    h_fileB = os.path.join(work_dir, "out.cen." + i + ".1D")

    if not os.path.exists(h_fileB):
        h_cmd = """
            3dToutcount -automask -fraction -polort {0} -legendre {1}+orig > {2}
            1deval -a {2} -expr "1-step(a-0.1)" > {3}
        """.format(
            pol,
            os.path.join(work_dir, i),
            h_fileA,
            h_fileB,
        )
        if test_mode:
            func_sbatch(h_cmd, 1, 1, 1, subj, sess, "outlier")
        else:
            func_afni(h_cmd)

    h_run = i.split("_")[0]
    for j in ["Forward", "Reverse"]:

        # create median datasets and masks
        h_blip_raw = os.path.join(work_dir, "blip_" + h_run + "_" + j)
        h_blip1 = os.path.join(work_dir, "tmp_blip_med_" + h_run + "_" + j)
        h_blip2 = os.path.join(work_dir, "tmp_blip_med_masked_" + h_run + "_" + j)

        if not os.path.exists(h_blip2 + "+orig.HEAD"):
            h_cmd = f"""
                3dTstat -median -prefix {h_blip1} {h_blip_raw}+orig
                3dAutomask -apply_prefix {h_blip2} {h_blip1}+orig
            """
            if test_mode:
                func_sbatch(h_cmd, 1, 1, 1, subj, sess, "fmap_med")
            else:
                func_afni(h_cmd)

    # comput midpoint warp, unwarp run data (solve for fall out), apply header
    h_q1 = os.path.join(work_dir, "tmp_blip_med_masked_" + h_run + "_Reverse")
    h_q2 = os.path.join(work_dir, "tmp_blip_med_masked_" + h_run + "_Forward")
    h_qout = os.path.join(work_dir, "tmp_blip_warp_" + h_run)

    h_blip_for = os.path.join(work_dir, "blip_" + h_run + "_Forward")
    h_run_epi = os.path.join(work_dir, i)

    if not os.path.exists(h_run_epi + "_blip+orig.HEAD"):
        h_cmd = f"""
            3dQwarp -plusminus -pmNAMES Rev For \
                -pblur 0.05 0.05 -blur -1 -1 \
                -noweight -minpatch 9 \
                -source {h_q1}+orig \
                -base {h_q2}+orig \
                -prefix {h_qout}

            3dNwarpApply -quintic -nwarp {h_qout}_For_WARP+orig \
                -source {h_run_epi}+orig -prefix {h_run_epi}_blip

            3drefit -atrcopy {h_blip_for}+orig IJK_TO_DICOM_REAL {h_run_epi}_blip+orig
        """
        if test_mode:
            func_sbatch(h_cmd, 1, 4, 2, subj, sess, "fmap_qwarp")
        else:
            func_afni(h_cmd)


# %%
# --- Step 3: Volume Registration
