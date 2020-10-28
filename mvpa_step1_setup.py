"""
Notes
"""
# %%
import subprocess
import fnmatch
import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from gp_step1_preproc import func_sbatch
from mvpa_step1_submit import task_dict as h_task_dict


# %%
def func_job(subj, subj_dir, decon_type, len_tr, task_dict, beh_dur, der_dir):
    """
    Step 1: Detrend

    1) Currently using "clean data" approach. Also, only using
        2GAM, not TENT data

        Should I just train on deconvolved sub-bricks?
    """

    # Work
    subj_num = subj.split("-")[1]
    mvpa_dir = os.path.join(der_dir, f"mvpa/sub{subj_num}")

    for phase in task_dict.keys():

        # get proper brick length
        #   REML appends an extra brick because
        #   "reasons". Account for AFNI's random
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
                            line.replace("#", "")
                            .split('"')[1]
                            .replace(" ", "")
                            .split(";")
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
            func_sbatch(h_cmd, 1, 4, 1, f"{subj_num}cle", subj_dir)

    # %%
    """
    Step 2: Organize MRI data

    1) Pymvpa expects certain files in certain locations
    """
    # pymvpa dirs
    py_dirs = ["BOLD", "anatomy", "model", "masks"]
    for i in py_dirs:
        if not os.path.exists(os.path.join(mvpa_dir, i)):
            os.makedirs(os.path.join(mvpa_dir, i))

    # anat
    if not os.path.exists(os.path.join(mvpa_dir, "anatomy/struct_ns.nii.gz")):
        h_cmd = f"module load afni-20.2.06 \n 3dcopy {subj_dir}/struct_ns+tlrc {mvpa_dir}/anatomy/struct_ns.nii.gz"
        h_job = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        print(h_job.communicate())

    # masks
    mask_dir = os.path.join(mvpa_dir, "masks", "orig")
    if not os.path.exists(mask_dir):
        os.makedirs(mask_dir)

    # pull, binarize priors
    atropos_dict = {2: "GMc", 4: "GMs"}
    atropos_dir = "/home/data/madlab/atlases/vold2_mni/priors_ACT"
    for i in atropos_dict:
        if not os.path.exists(
            os.path.join(subj_dir, f"tmp_{atropos_dict[i]}_bin.nii.gz")
        ):
            h_cmd = f"module load c3d/1.0.0 \n c3d {atropos_dir}/Prior{i}.nii.gz -thresh 0.3 1 1 0 -o {subj_dir}/tmp_{atropos_dict[i]}_bin.nii.gz"
            h_mask = subprocess.Popen(
                h_cmd, shell=True, stdout=subprocess.PIPE)
            out, err = h_mask.communicate()
            print(out, err)

    # make resampled GM-intersection mask
    if not os.path.exists(os.path.join(mask_dir, "GM_int_mask.nii.gz")):
        h_cmd = f"""
            module load c3d/1.0.0

            cd {subj_dir}
            c3d tmp_GMc_bin.nii.gz tmp_GMs_bin.nii.gz -add -o tmp_GM.nii.gz
            c3d tmp_GM.nii.gz -thresh 0.1 10 1 0 -o tmp_GM_bin.nii.gz

            3dfractionize -template CleanData_{phase}+tlrc -input tmp_GM_bin.nii.gz -prefix tmp_GM_res.nii.gz
            3dcalc -a tmp_GM_res.nii.gz -prefix tmp_GM_res_bin.nii.gz -expr 'step(a-3000)'

            3dcopy mask_epi_anat+tlrc tmp_mask_epi_anat.nii.gz
            c3d tmp_mask_epi_anat.nii.gz tmp_GM_res_bin.nii.gz -multiply -o {mask_dir}/GM_int_mask.nii.gz
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}msk", subj_dir)

    # BOLD - split into runs
    for count, phase in enumerate(task_dict.keys()):

        h_cmd = f"module load afni-20.2.06 \n 3dinfo -ntimes {subj_dir}/CleanData_{phase}+tlrc"
        h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        num_nvol = int(h_nvol.communicate()[0].decode("utf-8").strip())
        num_runs = len(scale_list)
        len_run = int(num_nvol / num_runs)

        beg_vol = 0
        end_vol = len_run - 1
        for run in range(1, num_runs + 1):
            bold_dir = os.path.join(mvpa_dir, f"BOLD/task00{count+1}_run00{run}")
            if not os.path.exists(bold_dir):
                os.makedirs(bold_dir)
            if not os.path.exists(os.path.join(bold_dir, "bold.nii.gz")):
                h_cmd = f"""
                    cd {subj_dir}
                    3dTcat -prefix tmp_run-{run}_{phase}_CleanData -tr {len_tr} "CleanData_{phase}+tlrc[{beg_vol}..{end_vol}]"
                    3dcopy tmp_run-{run}_{phase}_CleanData+tlrc {bold_dir}/bold.nii.gz
                """
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}spl", subj_dir)
            beg_vol += len_run
            end_vol += len_run

    # %%
    """
    Step 3: Organize timing files

    1) Based on timing files (e.g. tf_loc_face.txt) from
        gp_step2.R. Will generate needed attribute files.

        attributes.txt - one definition per volume
        cond00?.txt - each condition type has own onset time
            (onset duration amplitude(?) - 1.5 2 1)
            ref condtion_key.txt file for definitions

            Note - currently using float, problem?
    """
    for count, phase in enumerate(task_dict):
        cond_list = task_dict[phase]

        # split timing files into 1D
        len_sec = len_run * len_tr
        for cond in cond_list:
            h_cmd = f"module load afni-20.2.06 \n timing_tool.py -timing {subj_dir}/tf_{phase}_{cond}.txt -tr {len_tr} -stim_dur {beh_dur} -run_len {len_sec} -min_frac 0.3 -timing_to_1D {subj_dir}/tmp_tf_{phase}_{cond} -per_run_file"
            h_spl = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            print(h_spl.communicate())

        # make attribute files for e/run
        for run in range(1, num_runs + 1):

            # determine split tf files
            tf_list = [
                x
                for x in os.listdir(subj_dir)
                if fnmatch.fnmatch(x, f"tmp_tf_{phase}*_r0{run}.1D")
            ]
            tf_list.sort()

            # make df
            h_df = pd.concat(
                [
                    pd.read_csv(
                        os.path.join(subj_dir, item), names=[item.split("_")[3]]
                    )
                    for item in tf_list
                ],
                axis=1,
            )
            for cond in cond_list:
                h_df.loc[h_df[cond] > 0, "att"] = cond
            h_df = h_df.replace(np.nan, "base", regex=True)

            # subset df, write
            #   add col of 0s for reason?
            df_out = pd.DataFrame(h_df["att"])
            df_out["last"] = 0
            h_out = os.path.join(
                mvpa_dir, "BOLD", f"task00{count+1}_run00{run}", "attributes.txt"
            )
            np.savetxt(h_out, df_out.values, fmt="%s", delimiter=" ")

            # onset times - onset dur 1?
            model_dir = os.path.join(
                mvpa_dir, "model", "onsets", f"task00{count+1}_run00{run}"
            )
            if not os.path.exists(model_dir):
                os.makedirs(model_dir)
            for cc, cond in enumerate(cond_list):
                h_out = os.path.join(model_dir, f"cond00{cc+1}.txt")
                h_df = pd.read_csv(
                    os.path.join(subj_dir, f"tf_{phase}_{cond}.txt"),
                    sep="\t",
                    header=None,
                )
                df_out = pd.DataFrame(h_df.iloc[run - 1])
                df_out["dur"] = beh_dur
                df_out["last"] = 1
                df_out.drop(df_out.tail(1).index, inplace=True)
                np.savetxt(h_out, df_out.values, fmt="%s")


# receive arguments
def func_argparser():
    parser = ArgumentParser("Receive Bash args from wrapper")
    parser.add_argument("h_sub", help="Subject ID")
    parser.add_argument("h_dir", help="Subject Directory")
    parser.add_argument("h_dct", help="Decon Type")
    parser.add_argument("h_trl", help="TR Length")
    parser.add_argument("h_beh", help="Behavior Duration")
    parser.add_argument("h_der", help="Derivatives Directory")
    return parser


# %%
def main():

    # # For testing
    # h_subj = "sub-005"
    # h_subj_dir = "/scratch/madlab/nate_vCAT/derivatives/sub-005/ses-S1"
    # h_decon_type = "2GAM"
    # h_len_tr = 1.76
    # h_task_dict = {"loc": ["face", "scene", "num"]}
    # h_beh_dur = 1
    # h_der_dir = "/scratch/madlab/nate_vCAT/derivatives/sub-005"
    # func_job(h_subj, h_subj_dir, h_decon_type, h_len_tr, h_task_dict, h_beh_dur)

    args = func_argparser().parse_args()
    # print(args.h_sub, args.h_dir, args.h_dct, args.h_trl, h_task_dict, args.h_beh)
    func_job(
        args.h_sub, args.h_dir, args.h_dct, float(
            args.h_trl), h_task_dict, args.h_beh, args.h_der
    )


if __name__ == "__main__":
    main()

# %%
