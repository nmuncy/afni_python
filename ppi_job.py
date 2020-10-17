"""
Notes
"""

# %%
import os
import sys
import ast
import subprocess
import fnmatch
from step1_preproc import func_sbatch


# PPI decon function
def func_decon_ppi(run_files, mot_files, tf_dict, cen_file, h_str, h_type, ppi_dict):

    in_files = ""
    for fil in run_files:
        in_files += f"{fil} "

    reg_base = ""
    for c, mot in enumerate(mot_files):
        reg_base += f"-ortvec {mot} mot_demean_run{c+1} "

    c_beh = 0

    reg_beh = ""
    for beh in tf_dict:
        c_beh += 1
        if h_type == "dmBLOCK":
            reg_beh += f'-stim_times_AM1 {c_beh} {tf_dict[beh]} "dmBLOCK(1)" -stim_label {c_beh} {beh} '
        elif h_type == "GAM":
            reg_beh += (
                f'-stim_times {c_beh} {tf_dict[beh]} "GAM" -stim_label {c_beh} {beh} '
            )
        elif h_type == "2GAM":
            reg_beh += f'-stim_times {c_beh} {tf_dict[beh]} "TWOGAMpw(4,5,0.2,12,7)" -stim_label {c_beh} {beh} '

    for ts in ppi_dict:
        c_beh += 1
        reg_beh += f"-stim_file {c_beh} {ppi_dict[ts]} -stim_label {c_beh} {ts} "

    h_out = f"{h_str}_{h_type}_ppi"

    cmd_decon = f""" 3dDeconvolve \\
        -x1D_stop \\
        -input {in_files} \\
        -censor {cen_file} \\
        {reg_base} \\
        -polort A \\
        -float \\
        -local_times \\
        -num_stimts {c_beh} \\
        {reg_beh} \\
        -jobs 1 \\
        -x1D X.{h_out}.xmat.1D \\
        -xjpeg X.{h_out}.jpg \\
        -x1D_uncensored X.{h_out}.nocensor.xmat.1D \\
        -bucket {h_out}_stats \\
        -errts {h_out}_errts
    """
    return cmd_decon


# %%
def func_job(work_dir, subj, ses, phase, decon_type, seed_dict, stim_dur):
    # # for testing
    # work_dir = "/scratch/madlab/nate_vCAT/derivatives"
    # subj = "sub-006"
    # ses = "ses-S1"
    # phase = "vCAT"
    # decon_type = "2GAM"
    # seed_dict = {"LHC": "-24 -12 -22"}
    # stim_dur = 2

    # %%
    """
    Step 1: Clean Data

    Create "clean data" by removing effects of no interest
    (baseline regressors) from scaled data.
    """
    subj_dir = os.path.join(work_dir, subj, ses)

    # get TR
    h_cmd = (
        f"module load afni-20.2.06 \n 3dinfo -tr {subj_dir}/run-1_{phase}_scale+tlrc"
    )
    h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    len_tr = float(h_tr.communicate()[0].decode("utf-8").strip())

    # get proper brick length
    #   REML appends an extra brick because
    #   "reasons". Account for AFNIs random
    #   0-1 indexing
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -nv {subj_dir}/{phase}_{decon_type}_cbucket_REML+tlrc"
    h_len = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    len_wrong = h_len.communicate()[0].decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # list all scale files
    scale_list = [
        x.split(".")[0]
        for x in os.listdir(subj_dir)
        if fnmatch.fnmatch(x, f"*{phase}*scale+tlrc.HEAD")
    ]
    scale_list.sort()

    # make clean data
    if not os.path.exists(os.path.join(subj_dir, f"CleanData_{phase}+tlrc.HEAD")):

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
            3dTcat -prefix tmp_all_runs_{phase} -tr {len_tr} {" ".join(scale_list)}
            3dcalc -a tmp_all_runs_{phase}+tlrc -b tmp_effNoInt_{phase}+tlrc -expr 'a-b' -prefix CleanData_{phase}
        """
        func_sbatch(h_cmd, 1, 4, 1, "clean", subj_dir)

    # %%
    """
    Step 2: Seed Time Series

    Make HRF model, and seed from coordinates.
    Extract timeseries from ROI and upsample.
    Deconvolve HRF from timeseries (solve RHS).
    """
    # find smallest multiplier that returns int for resampling
    res_multiplier = 2
    status = False
    while not status:
        if ((len_tr * res_multiplier) % 2) == 0:
            status = True
        else:
            res_multiplier += 1

    # check for 1dUpsample
    if res_multiplier > 32:
        print(
            """
            Resample multiplier too high for use with 1dUpsample.
            Adjust code to continue.
            Exiting ...
        """
        )
        exit

    # make ideal HRF, use same model as deconvolution
    #   TENT not supported.
    if not os.path.exists(os.path.join(subj_dir, "HRF_model.1D")):

        if decon_type == "dmBLOCK":
            no_data = f"{14 + stim_dur} 1"
            hrf_model = "BLOCK(1,1)"
        elif decon_type == "GAM":
            no_data = "13 1"
            hrf_model = "GAM"
        elif decon_type == "2GAM":
            no_data = "19 1"
            hrf_model = "TWOGAMpw(4,5,0.2,12,7)"

        h_cmd = f"""
            cd {subj_dir}
            3dDeconvolve -polort -1 \
                -nodata {no_data} \
                -num_stimts 1 \
                -stim_times 1 1D:0 '{hrf_model}' \
                -x1D HRF_model.1D -x1D_stop
            1dUpsample {res_multiplier} HRF_model.1D > HRF_model_us.1D
        """
        func_sbatch(h_cmd, 1, 1, 1, "hrf", subj_dir)

    # get seed TS, solve RHS
    for key in seed_dict:

        # make seed, get TS
        if not os.path.exists(os.path.join(subj_dir, f"Seed_{key}+tlrc.HEAD")):
            h_cmd = f"""
                cd {subj_dir}
                echo {seed_dict[key]} | 3dUndump -xyz \
                    -srad 3 -master CleanData_{phase}+tlrc \
                    -prefix Seed_{key} -
                3dmaskave -quiet -mask Seed_{key}+tlrc CleanData_{phase}+tlrc > Seed_{key}_orig.1D
            """
            func_sbatch(h_cmd, 1, 1, 1, "mkseed", subj_dir)

        # solve RHS, then upsample
        #   I want to use -l2lasso, but I'm scared
        if not os.path.exists(os.path.join(subj_dir, f"Seed_{key}_neural_us.1D")):
            h_cmd = f"""
                cd {subj_dir}
                3dTfitter -RHS Seed_{key}_orig.1D -FALTUNG HRF_model.1D tmp.1D 012 0
                1dtranspose tmp.1D > Seed_{key}_neural.1D
                1dUpsample {res_multiplier} Seed_{key}_neural.1D > Seed_{key}_neural_us.1D
            """
            func_sbatch(h_cmd, 8, 1, 4, "faltung", subj_dir)

    # %%
    """
    Step 3: Make Behavior Time Series

    Make behavior vectors, extract behavior portions
        of seed neural timeseries.
    Convolve with HRF, then down sample
    """
    # list of timing files
    tf_list = [
        x for x in os.listdir(subj_dir) if fnmatch.fnmatch(x, f"tf_{phase}*.txt")
    ]

    # list of run length in seconds
    run_len = []
    for i in scale_list:
        h_cmd = f"module load afni-20.2.06 \n 3dinfo -ntimes {subj_dir}/{i}"
        h_vol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        num_vol = int(h_vol.communicate()[0].decode("utf-8").strip())
        run_len.append(num_vol * len_tr)

    for i in tf_list:
        h_beh = i.split(".")[0].split("_")[-1]

        # get upsampled behavior binary file
        if not os.path.exists(os.path.join(subj_dir, f"Beh_{h_beh}_bin.1D")):
            h_cmd = f"""
                cd {subj_dir}
                timing_tool.py -timing {i} -tr {len_tr} \
                    -stim_dur {stim_dur} -run_len {" ".join(map(str, run_len))} \
                    -min_frac 0.3 -timing_to_1D Beh_{h_beh}_bin.1D
                1dUpsample -1 {res_multiplier} Beh_{h_beh}_bin.1D > Beh_{h_beh}_us.1D
            """
            func_sbatch(h_cmd, 1, 1, 1, f"beh{h_beh}", subj_dir)

        # multiply beh bin file by clean neuro, convolve with HRF
        for key in seed_dict:
            if not os.path.exists(
                os.path.join(subj_dir, f"Seed_{key}_{h_beh}_bold_us.1D")
            ):
                h_cmd = f"""
                    cd {subj_dir}
                    1deval -a Seed_{key}_neural_us.1D -b Beh_{h_beh}_us.1D \
                        -expr 'a*b' > Seed_{key}_{h_beh}_neural_us.1D
                    waver -FILE {1 / res_multiplier} HRF_model_us.1D -input Seed_{key}_{h_beh}_neural_us.1D \
                        -numout {round((sum(run_len) / len_tr) * res_multiplier)} > Seed_{key}_{h_beh}_bold_us.1D
                """
                func_sbatch(h_cmd, 2, 1, 1, "behTS", subj_dir)

            # Downsample, bash > python
            if not os.path.exists(
                os.path.join(subj_dir, f"Final_{key}_{h_beh}_timeSeries.1D")
            ):
                h_cmd = f"""
                    cd {subj_dir}
                    cat Seed_{key}_{h_beh}_bold_us.1D | \
                        awk -v n={res_multiplier} 'NR%n==0' > Final_{key}_{h_beh}_timeSeries.1D
                """
                func_sbatch(h_cmd, 1, 1, 1, "dsts", subj_dir)

    # %%
    """
    Step 4: Rerun Deconvolution

    Same decon as before, add Seed and Behavior timeseries.
    """
    # Write decon script
    if not os.path.exists(
        os.path.join(subj_dir, f"X.{phase}_{decon_type}_ppi.xmat.1D")
    ):

        # determine motion list
        mot_list = [
            x
            for x in os.listdir(subj_dir)
            if fnmatch.fnmatch(x, f"mot_demean_{phase}.*.1D")
        ]
        mot_list.sort()

        # make timing file dict
        tf_dict = {}
        for i in tf_list:
            tmp = i.split("_")[-1]
            beh = tmp.split(".")[0]
            tf_dict[beh] = i

        # make ppi dict
        ppi_dict = {}
        for key in seed_dict:
            ppi_dict[key] = f"Seed_{key}_orig.1D"
            for beh in tf_dict:
                ppi_dict[f"{key}_{beh}"] = f"Final_{key}_{beh}_timeSeries.1D"

        # write decon script, generate matrices and REML_cmd
        decon_script = os.path.join(subj_dir, f"ppi_decon_{phase}_{decon_type}.sh")
        with open(decon_script, "w") as script:
            script.write(
                func_decon_ppi(
                    scale_list,
                    mot_list,
                    tf_dict,
                    f"censor_{phase}_combined.1D",
                    phase,
                    decon_type,
                    ppi_dict,
                )
            )

        # run decon script to generate matrices
        h_cmd = f"cd {subj_dir} \n source {decon_script}"
        func_sbatch(h_cmd, 1, 1, 1, "deconPPI", subj_dir)

    # run REML
    if not os.path.exists(
        os.path.join(subj_dir, f"{phase}_{decon_type}_ppi_stats_REML+tlrc.HEAD")
    ):
        h_cmd = f"cd {subj_dir} \n tcsh -x {phase}_{decon_type}_ppi_stats.REML_cmd -dsort {phase}_WMe_rall+tlrc"
        func_sbatch(h_cmd, 10, 4, 6, "reml", subj_dir)


def main():

    h_subj = str(sys.argv[1])
    h_ses = str(sys.argv[2])
    h_phase = str(sys.argv[3])
    h_decon_type = str(sys.argv[4])
    h_work_dir = str(sys.argv[5])
    h_seed_dict = ast.literal_eval(sys.argv[6])
    h_stim_dur = sys.argv[7]

    print(
        f"""
        subj = {h_subj}
        sess = {h_ses}
        phase = {h_phase}
        decon = {h_decon_type}
        work = {h_work_dir}
        seeds = {h_seed_dict}
        stim = {h_stim_dur}
    """
    )

    # func_job(h_work_dir, h_subj, h_ses, h_phase, h_decon_type, h_seed_dict, h_stim_dur)


if __name__ == "__main__":
    main()
