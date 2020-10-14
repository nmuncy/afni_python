"""
Notes
"""

# %%
import os
import subprocess
import fnmatch
from step1_preproc import func_sbatch

# %%
# for testing
work_dir = "/scratch/madlab/nate_vCAT/derivatives"
subj = "sub-006"
ses = "ses-S1"
phase = "vCAT"
decon_type = "2GAM"
seed_dict = {"LHC": "-24 -12 -22"}

# %%
"""
Step 1: clean data

Create "clean data" by removing effects of no interest
(baseline regressors) from scaled data.
"""
subj_dir = os.path.join(work_dir, subj, ses)

# get TR
h_cmd = f"module load afni-20.2.06 \n 3dinfo -tr {subj_dir}/run-1_{phase}_scale+tlrc"
h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
len_tr = h_tr.communicate()[0].decode("utf-8").strip()

# get proper brick length
#   REML appends an extra brick because
#   "reasons". Account for AFNIs random
#   0-1 indexing
h_cmd = f"module load afni-20.2.06 \n 3dinfo -nv {subj_dir}/{phase}_{decon_type}_cbucket_REML+tlrc"
h_len = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
len_wrong = h_len.communicate()[0].decode("utf-8").strip()
len_right = int(len_wrong) - 2

if not os.path.exists(os.path.join(subj_dir, f"CleanData_{phase}+tlrc.HEAD")):

    # list all scale files
    list_scale = [
        x.split(".")[0]
        for x in os.listdir(subj_dir)
        if fnmatch.fnmatch(x, f"*{phase}*scale+tlrc.HEAD")
    ]
    list_scale.sort()

    # list undesirable sub-bricks (those starting with Run or mot)
    no_int = []
    with open(os.path.join(subj_dir, f"X.{phase}_{decon_type}.xmat.1D")) as f:
        h_file = f.readlines()
        for line in h_file:
            if line.__contains__("ColumnLabels"):
                col_list = (
                    line.replace("#", "").split('"')[1].replace(" ", "").split(";")
                )
                for i, j in enumerate(col_list):
                    if fnmatch.fnmatch(j, "Run*") or fnmatch.fnmatch(j, "mot*"):
                        no_int.append(f"{str(i)}")

    # strip extra sub-brick, make clean data by removing
    #   effects of no interest from concatenated runs
    h_cmd = f"""
        cd {subj_dir}
        3dTcat -prefix tmp_{phase}_cbucket -tr {len_tr} "{phase}_{decon_type}_cbucket_REML+tlrc[0..{len_right}]"
        3dSynthesize -prefix tmp_effNoInt_{phase} -matrix X.{phase}_{decon_type}.xmat.1D \
            -cbucket tmp_{phase}_cbucket+tlrc -select {" ".join(no_int)} -cenfill nbhr
        3dTcat -prefix tmp_all_runs_{phase} -tr {len_tr} {" ".join(list_scale)}
        3dcalc -a tmp_all_runs_{phase}+tlrc -b tmp_effNoInt_{phase}+tlrc -expr 'a-b' -prefix CleanData_{phase}
    """
    func_sbatch(h_cmd, 1, 4, 1, "clean", subj_dir)

# %%
# make upsampled ideal 2GAM function
# h_mult = 2
# status = False
# while status:
#     if (float(len_tr) * h_mult) % 2 == 0:
#         status = True
#     else:
#         h_mult += 1
res_multiplier = 25

h_cmd = f"""
    3dDeconvolve -polort -1 \
        -nodata {int(float(len_tr) * res_multiplier)} 0.4 \
        -num_stimts 1 \
        -stim_times 1 1D:0 'TWOGAMpw(4,5,0.2,12,7)' \
        -x1D {subj_dir}/HRF_model.1D -x1D_stop
"""
func_sbatch(h_cmd, 1, 1, 1, "hrf", subj_dir)


# %%
# get seed TS, solve RHS
for key in seed_dict:

    # make seed, get TS and upsample
    if not os.path.exists(os.path.join(subj_dir, f"Seed_{key}+tlrc.HEAD")):
        h_cmd = f"""
            cd {subj_dir}
            echo {seed_dict[key]} | 3dUndump -xyz \
                -srad 3 -master CleanData_{phase}+tlrc \
                -prefix Seed_{key} -
            3dmaskave -quiet -mask Seed_{key}+tlrc \
                CleanData_{phase}+tlrc > Seed_{key}_ts_orig.1D
            1dUpsample {res_multiplier} Seed_{key}_ts_orig.1D > Seed_{key}_ts_us.1D
        """
        func_sbatch(h_cmd, 1, 1, 1, "mkseed", subj_dir)

    # solve RHS
    #   I want to use -l2lasso, but I'm scared
    if not os.path.exists(os.path.join(subj_dir, f"Seed_{key}_ts_neural.1D")):
        h_cmd = f"""
            cd {subj_dir}
            3dTfitter -RHS Seed_{key}_ts_us.1D \
                -FALTUNG HRF_model.1D tmp.1D 012 0
            1dtranspose tmp.1D > Seed_{key}_ts_neural.1D
        """
        func_sbatch(h_cmd, 2, 4, 1, "faltung", subj_dir)
