"""
Notes

Still only accepts one phase per session. Update this in the future.
"""

# %%
import os
import re
import sys
import subprocess
import fnmatch
from step1_preproc import func_sbatch


def func_pmBlock(tf_dict):
    h_reg_beh = ""
    for c_beh, beh in enumerate(tf_dict):
        h_reg_beh += f'-stim_times_AM1 {c_beh + 1} {tf_dict[beh]} "dmBLOCK(1)" -stim_label {c_beh + 1} {beh} '
    return h_reg_beh


def func_gam(tf_dict):
    h_reg_beh = ""
    for c_beh, beh in enumerate(tf_dict):
        h_reg_beh += f'-stim_times {c_beh + 1} {tf_dict[beh]} "GAM" -stim_label {c_beh + 1} {beh} '
    return h_reg_beh


def func_twoGam(df_dict):
    h_reg_beh = ""
    for c_beh, beh in enumerate(tf_dict):
        h_reg_beh += f'-stim_times {c_beh + 1} {tf_dict[beh]} "TWOGAMpw(4,5,0.2,12,7)" -stim_label {c_beh + 1} {beh} '
    return h_reg_beh


def func_decon(run_files, mot_files, tf_dict, cen_file, h_out, h_type):

    in_files = ""
    for fil in run_files:
        in_files += f"{fil.split('.')[0]} "

    reg_base = ""
    for c, mot in enumerate(mot_files):
        reg_base += f"-ortvec {mot} mot_demean_run{c+1} "

    if h_type == "pmBLOCK":
        reg_beh = func_pmBlock(tf_dict)
    elif h_type == "GAM":
        reg_beh = func_gam(tf_dict)
    elif h_type == "2GAM":
        reg_beh = func_twoGam(tf_dict)

    cmd_decon = f""" 3dDeconvolve \\
        -x1D_stop \\
        -input {in_files} \\
        -censor {cen_file} \\
        {reg_base} \\
        -polort A \\
        -float \\
        -local_times \\
        -num_stimts {len(tf_dict.keys())} \\
        {reg_beh} \\
        -jobs 1 \\
        -x1D X.{h_out}.xmat.1D \\
        -xjpeg X.{h_out}.jpg \\
        -x1D_uncensored X.{h_out}.nocensor.xmat.1D \\
        -bucket {h_out}_stats -errts {h_out}_errts
    """
    return cmd_decon


# For testing
subj = "sub-006"
sess = "ses-S1"
phase = "vCAT"

par_dir = "/scratch/madlab/nate_vCAT"
data_dir = os.path.join(par_dir, "dset", subj, sess)
work_dir = os.path.join(par_dir, "derivatives", subj, sess)


# %%
"""
Step 1: Make motion regressors
"""
run_list = [
    x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, f"*{phase}_scale+tlrc.HEAD")
]
num_run = len(run_list)
run_list.sort()

h_cmd = f"""
    cd {work_dir}

    cat dfile.run-*_{phase}.1D > dfile_rall_{phase}.1D

    # files: de-meaned, motion params (per phase)
    1d_tool.py -infile dfile_rall_{phase}.1D \
        -set_nruns {num_run} \
        -demean \
        -write motion_demean_{phase}.1D

    1d_tool.py -infile dfile_rall_{phase}.1D \
        -set_nruns {num_run} \
        -derivative \
        -demean \
        -write motion_deriv_{phase}.1D

    1d_tool.py -infile motion_demean_{phase}.1D \
        -set_nruns {num_run} \
        -split_into_pad_runs mot_demean_{phase}

    1d_tool.py -infile dfile_rall_{phase}.1D \
        -set_nruns {num_run} \
        -show_censor_count \
        -censor_prev_TR \
        -censor_motion 0.3 motion_{phase}

    # determine censor
    cat out.cen.run-*{phase}.1D > outcount_censor_{phase}.1D
    1deval -a motion_{phase}_censor.1D -b outcount_censor_{phase}.1D \
        -expr "a*b" > censor_{phase}_combined.1D
"""
if not os.path.exists(os.path.join(work_dir, f"censor_{phase}_combined.1D")):
    func_sbatch(h_cmd, 1, 1, 1, "motion", work_dir)


# %%
"""
Step 2: Deconvolve
"""

# Get motion files
mot_list = [
    x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, f"mot_demean_{phase}.*.1D")
]
mot_list.sort()


# Get timing files - all are named *_phase_beh.txt
#   e.g. tf_vCAT_hit.txt
tf_list = [x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, f"*{phase}*.txt")]

# make timing file dictionary
tf_dict = {}
for i in tf_list:
    tmp = i.split("_")[-1]
    beh = tmp.split(".")[0]
    tf_dict[beh] = i

# write decon script, generate matrices and REML_cmd
decon_script = os.path.join(work_dir, f"decon_{phase}.sh")
with open(decon_script, "w") as script:
    script.write(
        func_decon(
            run_list, mot_list, tf_dict, f"censor_{phase}_combined.1D", phase, "pmBLOCK"
        )
    )

if not os.path.exists(os.path.join(work_dir, f"X.{phase}.xmat.1D")):
    h_cmd = f"cd {work_dir} \n source {decon_script}"
    func_sbatch(h_cmd, 1, 1, 1, "decon", work_dir)

# generate WM timeseries
if not os.path.exists(os.path.join(work_dir, f"{phase}_WMe_rall+tlrc.HEAD")):
    h_cmd = f"""
        cd {work_dir}
        3dTcat -prefix tmp_allRuns_{phase} run-*{phase}_scale+tlrc.HEAD
        3dcalc -a tmp_allRuns_{phase}+tlrc -b final_mask_WM_eroded+tlrc \
            -expr 'a*bool(b)' -datum float -prefix tmp_allRuns_{phase}_WMe
        3dmerge -1blur_fwhm 20 -doall -prefix {phase}_WMe_rall tmp_allRuns_{phase}_WMe+tlrc
    """
    func_sbatch(h_cmd, 1, 1, 4, "wmts", work_dir)

# run REML
if not os.path.exists(os.path.join(work_dir, f"{phase}_stats_REML+tlrc.HEAD")):
    h_cmd = (
        f"cd {work_dir} \n tcsh -x {phase}_stats.REML_cmd -dsort {phase}_WMe_rall+tlrc"
    )
    func_sbatch(h_cmd, 10, 6, 4, "reml", work_dir)
