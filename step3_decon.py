"""
Notes

Timing files should be named "tf_phase_behavior.txt"
    e.g. tf_vCAT_Hit.txt, or tf_test_FA.txt

If using dmBLOCK option, duration should be married
    to start time.

GAM, 2GAM base models not using duration - just use
    BLOCK!

TODO:
    1) Finish writing/testing for multiple phases
    2) Test different basis models
"""

# %%
import os
import re
import sys
import subprocess
import fnmatch
from step1_preproc import func_sbatch


def func_decon(run_files, mot_files, tf_dict, cen_file, h_str, h_type):

    in_files = ""
    for fil in run_files:
        in_files += f"{fil.split('.')[0]} "

    reg_base = ""
    for c, mot in enumerate(mot_files):
        reg_base += f"-ortvec {mot} mot_demean_run{c+1} "

    reg_beh = ""
    for c_beh, beh in enumerate(tf_dict):
        if h_type == "dmBLOCK":
            reg_beh += f'-stim_times_AM1 {c_beh + 1} {tf_dict[beh]} "dmBLOCK(1)" -stim_label {c_beh + 1} {beh} '
        elif h_type == "GAM":
            reg_beh += f'-stim_times {c_beh + 1} {tf_dict[beh]} "GAM" -stim_label {c_beh + 1} {beh} '
        elif h_type == "2GAM":
            reg_beh += f'-stim_times {c_beh + 1} {tf_dict[beh]} "TWOGAMpw(4,5,0.2,12,7)" -stim_label {c_beh + 1} {beh} '

    h_out = f"{h_str}_{h_type}"

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


# %%
def func_job(phase, decon_type, work_dir):

    # # For testing
    # subj = "sub-006"
    # sess = "ses-S1"
    # phase = "vCAT"
    # decon_type = "GAM"

    # par_dir = "/scratch/madlab/nate_vCAT"
    # work_dir = os.path.join(par_dir, "derivatives", subj, sess)

    """
    Step 1: Make motion regressors

    Creates motion, demean, and derivative files. Demeaned
        are the ones used in the deconvolution.

    Censor file is combination of 2 things:
        1) Censors based on >0.3 rot/translation relative to previous
            volume. Previous volume censored also.
        2) Which volumes had >10% outlier voxels.
    """
    run_list = [
        x
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, f"*{phase}_scale+tlrc.HEAD")
    ]
    num_run = len(run_list)
    run_list.sort()

    h_cmd = f"""
        cd {work_dir}

        cat dfile.run-*_{phase}.1D > dfile_rall_{phase}.1D

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

        cat out.cen.run-*{phase}.1D > outcount_censor_{phase}.1D
        1deval -a motion_{phase}_censor.1D -b outcount_censor_{phase}.1D \
            -expr "a*b" > censor_{phase}_combined.1D
    """
    if not os.path.exists(os.path.join(work_dir, f"censor_{phase}_combined.1D")):
        func_sbatch(h_cmd, 1, 1, 1, "motion", work_dir)

    # %%
    """
    Step 2: Deconvolve

    Uses 3dDeconvolve to generate matrix. 3dREMLfit is then used
        to do a GLS with an ARMA function.

    White matter time series is used as a nuissance regressor.

    Deconvolve script written for review.

    Base models include pmBLOCK, GAM, and TWOGAMpw.
    """

    # Get motion files
    mot_list = [
        x
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, f"mot_demean_{phase}.*.1D")
    ]
    mot_list.sort()

    # Get timing files - all are named tf_phase_beh.txt
    #   e.g. tf_vCAT_hit.txt
    tf_list = [
        x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, f"tf*{phase}*.txt")
    ]

    # make timing file dictionary
    tf_dict = {}
    for i in tf_list:
        tmp = i.split("_")[-1]
        beh = tmp.split(".")[0]
        tf_dict[beh] = i

    # write decon script, generate matrices and REML_cmd
    decon_script = os.path.join(work_dir, f"decon_{phase}_{decon_type}.sh")
    with open(decon_script, "w") as script:
        script.write(
            func_decon(
                run_list,
                mot_list,
                tf_dict,
                f"censor_{phase}_combined.1D",
                phase,
                decon_type,
            )
        )

    # run decon script to generate matrices
    if not os.path.exists(os.path.join(work_dir, f"X.{phase}_{decon_type}.xmat.1D")):
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
    if not os.path.exists(
        os.path.join(work_dir, f"{phase}_{decon_type}_stats_REML+tlrc.HEAD")
    ):
        h_cmd = f"cd {work_dir} \n tcsh -x {phase}_{decon_type}_stats.REML_cmd -dsort {phase}_WMe_rall+tlrc"
        func_sbatch(h_cmd, 10, 6, 4, "reml", work_dir)


def main():

    subj = str(sys.argv[1])
    sess = str(sys.argv[2])
    h_phase = str(sys.argv[3])
    h_decon_type = str(sys.argv[4])
    h_deriv_dir = str(sys.argv[5])

    h_work_dir = os.path.join(h_deriv_dir, subj, sess)
    func_job(h_phase, h_decon_type, h_work_dir)


if __name__ == "__main__":
    main()
