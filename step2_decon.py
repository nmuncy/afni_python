"""
Notes

Still only accepts one phase per session. Update this in the future.
"""

# %%
import os
import sys
import subprocess
import fnmatch
from step1_preproc import func_sbatch


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


def func_pmBlock(tf_dict):
    reg_beh = ""
    for c_beh, beh in enumerate(tf_dict):
        reg_beh += f"""
            -stim_times_AM2 {c_beh + 1} {tf_dict[beh]} \"dmBLOCK(1)\" \
            -stim_label {c_beh + 1} {beh}
        """
    return reg_beh


# tf_dict = {"BehA": "Timing_fileA.txt"}
def func_decon(run_files, mot_files, tf_dict, cen_file, h_out, h_type):

    # make reg_base
    reg_base = ""
    for c, mot in enumerate(mot_files):
        reg_base += f"-ortvec {mot} mot_demean_run{c+1} "

    # make reg_beh
    if h_type == "pmBlock":
        reg_beh = func_pmBlock(tf_dict)

    cmd_decon = f"""
        3dDeconvolve \
        -x1D_stop \
        -input {run_files} \
        -censor {cen_file} \
        {reg_base} \
        -polort A \
        -float \
        -num_stimts {len(tf_dict.keys())} \
        {reg_beh} \
        -jobs 1 \
        -x1D X.{h_out}.xmat.1D \
        -xjpeg X.{h_out}.jpg \
        -x1D_uncensored X.{h_out}.nocensor.xmat.1D \
        -bucket {h_out}_stats -errts {h_out}_errts
    """
    return cmd_decon


# Get motion files
mot_list = [
    x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, f"mot_demean_{phase}.*.1D")
]
mot_list.sort()


# # Get timing files
# tf_list = [
#     x
#     for x in os.listdir(os.path.join(work_dir, "timing_files"))
#     if fnmatch.fnmatch(x, f"*{phase}*.txt")
# ]

# make timing file dictionary
tf_dict = {}


# write decon script
decon_script = os.path.join(work_dir, f"decon_{phase}.sh")
with open(decon_script, "w") as script:
    script.write(
        func_decon(
            run_list, mot_list, tf_dict, f"censor_{phase}_combined.1D", phase, "pmBLOCK"
        )
    )

h_cmd = f"cd {work_dir} \n source decon_{phase}.sh"
func_sbatch(h_cmd, 1, 1, 1, "decon", work_dir)
