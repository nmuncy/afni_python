# %%
# --- Notes
#
# 1) only written to accept one experiment phase (e.g. Study, Test) per session, assumes
#       all epi scans in session pertain to phase
#
# TODO:
#   1) Update wait function to wait for job name rather than job number
#   2) Update template?
#   3) Test various functions/syntax


import json
import os
import sys
import subprocess
import fnmatch
import math
import time
import codecs


# %%
# --- Step 0: Set up
#
# Set/receive arguments, set variables for data, working, and atlas directories,
#   make working dir if needed, and set up general functions.

test_mode = True

if test_mode:
    subj = "sub-005"
    sess = "ses-S1"
    phase = "vCAT"
else:
    subj = str(sys.argv[1])
    sess = str(sys.argv[2])
    phase = str(sys.argv[3])

par_dir = "/scratch/madlab/nate_vCAT"
data_dir = os.path.join(par_dir, "dset", subj, sess)
work_dir = os.path.join(par_dir, "derivatives", subj, sess)
atlas_dir = "TODO"

if not os.path.exists(work_dir):
    os.makedirs(work_dir)


# Submit jobs to slurm
#   Note: len(h_str) < 8
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_sub, h_ses, h_str):

    full_name = h_sub + "_" + h_ses + "_" + h_str
    sbatch_job = f"""
        sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load afni-20.2.06 \n {command}"
    """

    print(sbatch_job)
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id, error = sbatch_response.communicate()

    while_count = 0
    status = False
    while not status:

        check_cmd = "squeue -u $(whoami)"
        sq_check = subprocess.Popen(check_cmd, shell=True, stdout=subprocess.PIPE)
        out_lines, err_lines = sq_check.communicate()
        b_decode = out_lines.decode("utf-8")

        if h_str not in b_decode:
            status = True

        if not status:
            while_count += 1
            print(f"Wait count for sbatch job {h_str}: ", while_count)
            time.sleep(3)
    print(f'Sbatch job "{h_str}" finished')


# Submit afni subprocess
# TODO: This function is untested
def func_afni(h_cmd):
    cmd = f"module load afni-20.2.06 \n {h_cmd}"
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
#
# Make epi_dict for later use, like pulling header info
#   or looping through runs.

# struct
struct_nii = os.path.join(data_dir, "anat", "{}_{}_T1w.nii.gz".format(subj, sess))
struct_raw = os.path.join(work_dir, "struct+orig")
if not os.path.exists(struct_raw + ".HEAD"):
    h_cmd = "3dcopy {} {}".format(struct_nii, struct_raw)
    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "struct")
    else:
        func_afni(h_cmd)


# epi - only task, not rsFMRI
epi_list = [
    epi
    for epi in os.listdir(os.path.join(data_dir, "func"))
    if fnmatch.fnmatch(epi, "*task-vCAT*.nii.gz")
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

# %%
# fmap
#   In EMU, account for potential multiple fmaps, use IntendedFor
json_list = [
    x
    for x in os.listdir(os.path.join(data_dir, "fmap"))
    if fnmatch.fnmatch(x, "*.json")
]

fmap_list = []
for i in json_list:
    with open(os.path.join(data_dir, "fmap", i)) as j:
        h_json = json.load(j)
        if "IntendedFor" not in h_json:
            fmap_list.append(i.split(".")[0] + ".nii.gz")
        else:
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


# Make fmap for each epi run
epi_len = len(epi_list) + 1
if len(fmap_list) == 2:
    for i in range(1, epi_len):
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

for i in epi_dict.keys():

    # determine polort arg, file strings, find outliers
    h_fileA = os.path.join(work_dir, "outcount." + i + ".1D")
    h_fileB = os.path.join(work_dir, "out.cen." + i + ".1D")

    if not os.path.exists(h_fileB):
        h_cmd = """
            len_tr=`3dinfo -tr {0}+orig`
            pol_time=$(echo $(echo $hold*$len_tr | bc)/150 | bc -l)
            pol=$((1 + `printf "%.0f" $pol_time`))
            3dToutcount -automask -fraction -polort $pol -legendre {0}+orig > {1}
            1deval -a {1} -expr "1-step(a-0.1)" > {2}
        """.format(
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
            func_sbatch(h_cmd, 1, 4, 2, subj, sess, "qwarp")
        else:
            func_afni(h_cmd)


# %%
# --- Step 3: Make volreg base
#
# 1) The volreg base (epi_vr_base) is the single volume in entire experiment
#   phase with the smallest number of outlier volumes.
#
# Note: If data is not time shifted, do so before determining volreg_base.
#
# Also, it was easier to write a small bash script due to the multiple
#   afni commands.


out_list = [
    os.path.join(work_dir, x)
    for x in os.listdir(work_dir)
    if fnmatch.fnmatch(x, "outcount.*.1D")
]
out_all = os.path.join(work_dir, "outcount_all.1D")
with open(out_all, "w") as outfile:
    for i in out_list:
        with open(i) as infile:
            outfile.write(infile.read())

if not os.path.exists(os.path.join(work_dir, "epi_vr_base+orig.HEAD")):
    vr_script = os.path.join(work_dir, "do_volregBase.sh")
    with open(vr_script, "w") as script:
        script.write(
            """
            #!/bin/bash
            cd {}

            unset tr_counts block
            numRuns=0; for i in run*_vCAT+orig.HEAD; do
                hold=`3dinfo -ntimes ${{i%.*}}`
                tr_counts+="$hold "
                block[$numRuns]=${{i%+*}}
                let numRuns=$[$numRuns+1]
            done
            # echo ${{block[@]}}

            minindex=`3dTstat -argmin -prefix - outcount_all.1D\\'`
            ovals=(`1d_tool.py -set_run_lengths $tr_counts -index_to_run_tr $minindex`)
            minoutrun=${{ovals[0]}}
            minouttr=${{ovals[1]}}
            # echo $minoutrun $minouttr

            c=0; for ((d=1; d <= $numRuns; d++)); do
                if [ 0$d == $minoutrun ]; then
                    baseRun=${{block[$c]}}
                fi
                let c=$[$c+1]
            done
            # echo $baseRun

            3dbucket -prefix epi_vr_base ${{baseRun}}_blip+orig"[${{minouttr}}]"
            """.format(
                work_dir
            )
        )

    h_cmd = f"source {vr_script}"
    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "vrbase")
    else:
        func_afni(h_cmd)

# %%
# --- Step 4: Calc, Perfrom normalization
#
# This step will perform the rigid alignments of T1-EPI (A)
#   EPI-EPI base volume (B), and non-linear diffeomorphic of T1-Template (C).
#
# It will then concatenate these warp matrices, and warp EPI data from
#   raw/native space to template space via W=A'+B+C. Thus, only one
#   interpolation of the epi data occurs.

# Calculate T1-EPI rigid, T1-Template diffeo
h_cmd = f"""
    template=~/bin/Templates/vold2_mni/vold2_mni_brain+tlrc
    cd {work_dir}

    align_epi_anat.py \\
        -anat2epi \\
        -anat struct+orig \\
        -save_skullstrip \\
        -suffix _al_junk \\
        -epi epi_vr_base+orig \\
        -epi_base 0 \\
        -epi_strip 3dAutomask \\
        -cost lpc+ZZ \\
        -volreg off \\
        -tshift off

    auto_warp.py -base $template -input struct_ns+orig -skull_strip_input no
    3dbucket -prefix struct_ns awpy/struct_ns.aw.nii*
    cp awpy/anat.un.aff.Xat.1D .
    cp awpy/anat.un.aff.qw_WARP.nii .
"""
if not os.path.exists(os.path.join(work_dir, "struct_ns+tlrc.HEAD")):
    if test_mode:
        func_sbatch(h_cmd, 4, 4, 4, subj, sess, "diffeo")
    else:
        func_afni(h_cmd)

# Calculate volreg for e/run
for i in epi_dict.keys():

    h_cmd = f"""
        cd {work_dir}

        3dvolreg -verbose \
        -zpad 1 \
        -base epi_vr_base+orig \
        -1Dfile dfile.{i}.1D \
        -prefix {i}_volreg \
        -cubic \
        -1Dmatrix_save mat.{i}.vr.aff12.1D \
        {i}_blip+orig
    """
    if not os.path.exists(os.path.join(work_dir, f"mat.{i}.vr.aff12.1D")):
        if test_mode:
            func_sbatch(h_cmd, 1, 1, 1, subj, sess, "volreg")
        else:
            func_afni(h_cmd)


# %%
# Concat calcs, warp EPI. Make mask
for i in epi_dict.keys():

    h_cmd = f"""
        cd {work_dir}

        cat_matvec -ONELINE \
        anat.un.aff.Xat.1D \
        struct_al_junk_mat.aff12.1D -I \
        mat.{i}.vr.aff12.1D > mat.{i}.warp.aff12.1D

        gridSize=`3dinfo -di {i}+orig`

        3dNwarpApply -master struct_ns+tlrc \
        -dxyz $gridSize \
        -source {i}_blip+orig \
        -nwarp "anat.un.aff.qw_WARP.nii mat.{i}.warp.aff12.1D" \
        -prefix {i}_warp

        3dcalc -overwrite -a {i}_blip+orig -expr 1 -prefix tmp_{i}_mask

        3dNwarpApply -master struct_ns+tlrc \
        -dxyz $gridSize \
        -source tmp_{i}_mask+orig \
        -nwarp "anat.un.aff.qw_WARP.nii mat.{i}.warp.aff12.1D" \
        -interp cubic \
        -ainterp NN -quiet \
        -prefix {i}_mask_warped

        3dTstat -min -prefix tmp_{i}_min {i}_mask_warped+tlrc
    """
    if not os.path.exists(os.path.join(work_dir, f"{i}_warp+tlrc.HEAD")):
        if test_mode:
            func_sbatch(h_cmd, 1, 4, 4, subj, sess, "warp")
        else:
            func_afni(h_cmd)

# %%
# Determine minimum value, make mask
h_cmd = f"""
    cd {work_dir}
    3dMean -datum short -prefix tmp_mean tmp_run-{{1..{len(epi_dict.keys())}}}_{phase}_min+tlrc
    3dcalc -a tmp_mean+tlrc -expr 'step(a-0.999)' -prefix {phase}_minVal_mask
"""
# print(h_cmd)
if not os.path.exists(os.path.join(work_dir, f"{phase}_minVal_mask+tlrc.HEAD")):
    if test_mode:
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "minVal")
    else:
        func_afni(h_cmd)


# %%
# --- Step 5: Scale data
