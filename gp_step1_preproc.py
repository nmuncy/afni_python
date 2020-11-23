# %%
"""
Notes

1) Script will generate named sbatch subprocesses named for each sub-step.
    Output for each sbatch step is written to work_dir and prepeneded
    with "sbatch_writeOut"

2) Pre-processing steps include copying data to derivatives, distortion correcting,
    volume registration, normalization, generating intersection masks, blurring, and
    scaling the data.

3) Written in Python 3.8, has afni and c3d dependencies.
"""

import json
import os
import subprocess
import fnmatch
import math
from argparse import ArgumentParser
from gp_step0_dcm2nii import func_sbatch


# make a list of all EPI scans of phase X
def func_epi_list(phase, h_dir):
    h_list = [
        x.split("+")[0]
        for x in os.listdir(h_dir)
        if fnmatch.fnmatch(x, f"run-*_{phase}+orig.HEAD")
    ]
    h_list.sort()
    return h_list


# receive arguments
def func_argparser():
    parser = ArgumentParser("Receive Bash args from wrapper")
    parser.add_argument("h_sub", help="Subject ID")
    parser.add_argument("h_ses", help="Session")
    parser.add_argument("h_blt", type=int, help="Blip Toggle")
    parser.add_argument("h_par", help="Parent Directory")
    parser.add_argument("h_phl", nargs="+", help="Phase List")
    return parser


# %%
def func_preproc(data_dir, work_dir, subj, sess, phase_list, blip_tog):

    """
    Step 1: Copy data into work_dir

    1) Get func, anat, fmap data. Rename appropriately.

    2) To account for different num of fmaps, will
        produce AP, PA fmap per run.
    """

    # # For testing
    # subj = "sub-005"
    # sess = "ses-S1"
    # phase_list = ["loc"]
    # blip_tog = 1

    # par_dir = "/scratch/madlab/nate_vCAT"
    # data_dir = os.path.join(par_dir, "dset", subj, sess)
    # work_dir = os.path.join(par_dir, "derivatives", subj, sess)

    # if not os.path.exists(work_dir):
    #     os.makedirs(work_dir)

    # Start
    subj_num = subj.split("-")[1]

    # struct
    struct_nii = os.path.join(data_dir, "anat", f"{subj}_{sess}_T1w.nii.gz")
    struct_raw = os.path.join(work_dir, "struct+orig")
    if not os.path.exists(struct_raw + ".HEAD"):
        h_cmd = f"3dcopy {struct_nii} {struct_raw}"
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}t1w", work_dir)

    # epi - keep runs sorted by phase
    for phase in phase_list:
        epi_list = [
            epi
            for epi in os.listdir(os.path.join(data_dir, "func"))
            if fnmatch.fnmatch(epi, f"*task-{phase}*.nii.gz")
        ]
        epi_list.sort()
        for i in range(len(epi_list)):
            epi_raw = os.path.join(work_dir, f"run-{1 + i}_{phase}+orig")
            epi_nii = os.path.join(data_dir, "func", epi_list[i])
            if not os.path.exists(epi_raw + ".HEAD"):
                h_cmd = f"3dcopy {epi_nii} {epi_raw}"
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}epi", work_dir)

    # %%
    # fmap
    if blip_tog == 1:

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
                    for phase in phase_list:
                        epi_list = [
                            epi
                            for epi in os.listdir(os.path.join(data_dir, "func"))
                            if fnmatch.fnmatch(epi, f"*task-{phase}*.nii.gz")
                        ]
                        for k in epi_list:
                            h_epi = os.path.join(sess, "func", k)
                            if h_epi in h_json["IntendedFor"]:
                                h_fmap = i.split(".")[0] + ".nii.gz"
                                if h_fmap not in fmap_list:
                                    fmap_list.append(h_fmap)

        # copy fmap function
        def func_fmap(h_file, h_run):

            fmap_nii = os.path.join(data_dir, "fmap", h_file)
            h_dir = h_file.split("-")[4].lstrip().split("_")[0]
            enc_dir = "Forward" if h_dir == "AP" else "Reverse"
            fmap_raw = os.path.join(work_dir, f"blip_run-{h_run}_{enc_dir}+orig")

            if not os.path.exists(fmap_raw + ".HEAD"):
                h_cmd = f"3dcopy {fmap_nii} {fmap_raw}"
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}fmap", work_dir)

        # Make fmap for each epi run
        #   TODO: Test this more (if len(fmap_list) != 2)
        for phase in phase_list:
            epi_list = [
                epi
                for epi in os.listdir(os.path.join(data_dir, "func"))
                if fnmatch.fnmatch(epi, f"*task-{phase}*.nii.gz")
            ]
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
    """
    Step 2: Detect outliers voxels, blip correct

    1) Determine the proportion of voxels per volume that have outlier signal.
        Censor volumes that exceed limit.

    2) Correct for signal fallout using fmap. This approach is taken from
        afni_proc. It uses the fmap to "unwarp" the run epi.
    """

    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        for run in epi_list:

            # calc tr length
            h_cmd = (
                f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -tr {run}+orig"
            )
            h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_len_tr = h_tr.communicate()[0]
            len_tr = h_len_tr.decode("utf-8").strip()

            # calc number of volumes
            h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -ntimes {run}+orig"
            h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_nvol_out = h_nvol.communicate()[0]
            num_tr = h_nvol_out.decode("utf-8").strip()

            # determine polort argument
            pol = 1 + math.ceil((float(len_tr) * float(num_tr)) / 150)

            # determine percentage outliers for e/volume
            if not os.path.exists(os.path.join(work_dir, f"out.cen.{run}.1D")):
                h_cmd = f"""
                    cd {work_dir}
                    3dToutcount -automask -fraction -polort {pol} \
                        -legendre {run}+orig > outcount.{run}.1D
                    1deval -a outcount.{run}.1D -expr '1-step(a-0.05)' > out.cen.{run}.1D
                """
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}out", work_dir)

            if blip_tog == 1:
                h_run = run.split("_")[0]
                for j in ["Forward", "Reverse"]:

                    # create median datasets and masks
                    if not os.path.exists(
                        os.path.join(
                            work_dir, f"tmp_blip_med_masked_{h_run}_{j}+orig.HEAD"
                        )
                    ):
                        h_cmd = f"""
                            cd {work_dir}
                            3dTstat -median -prefix tmp_blip_med_{h_run}_{j} blip_{h_run}_{j}+orig
                            3dAutomask -apply_prefix tmp_blip_med_masked_{h_run}_{j} \
                                tmp_blip_med_{h_run}_{j}+orig
                        """
                        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}med", work_dir)

                # comput midpoint warp, unwarp run data (solve for fall out), apply header
                if not os.path.exists(os.path.join(work_dir, f"{run}_blip+orig.HEAD")):
                    h_cmd = f"""
                        cd {work_dir}

                        3dQwarp -plusminus -pmNAMES Rev For \
                            -pblur 0.05 0.05 -blur -1 -1 \
                            -noweight -minpatch 9 \
                            -source tmp_blip_med_masked_{h_run}_Reverse+orig \
                            -base tmp_blip_med_masked_{h_run}_Forward+orig \
                            -prefix blip_{h_run}

                        3dNwarpApply -quintic -nwarp blip_{h_run}_For_WARP+orig \
                            -source {run}+orig -prefix {run}_blip

                        3drefit -atrcopy blip_{h_run}_Forward+orig IJK_TO_DICOM_REAL {run}_blip+orig
                    """
                    func_sbatch(h_cmd, 1, 2, 4, f"{subj_num}qwa", work_dir)

    # %%
    """
    Step 3: Make volreg base

    1) The volreg base (epi_vrBase) is the single volume in entire experiment
        with the smallest number of outlier volumes.

    Note: If data is not time shifted, do so before determining volreg_base.
        Also, it was easier to write a small bash script due to the multiple
        afni commands.
    """

    out_list = [
        os.path.join(work_dir, x)
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, "outcount.*.1D")
    ]
    out_list.sort()
    out_all = os.path.join(work_dir, "outcount_all.1D")
    with open(out_all, "w") as outfile:
        for i in out_list:
            with open(i) as infile:
                outfile.write(infile.read())

    if not os.path.exists(os.path.join(work_dir, "epi_vrBase+orig.HEAD")):
        vr_script = os.path.join(work_dir, "do_volregBase.sh")
        if blip_tog == 1:
            run_str = "blip+orig"
        else:
            run_str = "+orig"
        with open(vr_script, "w") as script:
            script.write(
                """
                #!/bin/bash
                cd {}
                unset tr_counts block
                numRuns=0; for i in run-*_{}.HEAD; do
                    hold=`3dinfo -ntimes ${{i%.*}}`
                    tr_counts+="$hold "
                    block[$numRuns]=${{i%+*}}
                    let numRuns=$[$numRuns+1]
                done

                minindex=`3dTstat -argmin -prefix - outcount_all.1D\\'`
                ovals=(`1d_tool.py -set_run_lengths $tr_counts -index_to_run_tr $minindex`)
                minoutrun=${{ovals[0]}}
                minouttr=${{ovals[1]}}

                c=0; for ((d=1; d <= $numRuns; d++)); do
                    if [ 0$d == $minoutrun ]; then
                        baseRun=${{block[$c]}}
                    fi
                    let c=$[$c+1]
                done

                3dbucket -prefix epi_vrBase ${{baseRun}}+orig"[${{minouttr}}]"
                """.format(
                    work_dir, run_str
                )
            )

        h_cmd = f"source {vr_script}"
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}vrb", work_dir)

    # %%
    """
    Step 4: Calc, Perfrom normalization

    1) This step will perform the rigid alignments of T1-EPI (A)
        EPI-EPI base volume (B), and non-linear diffeomorphic of
        T1-Template (C). Note blip distortion map (D).

    2) It will then concatenate these warp matrices, and warp EPI data from
        raw/native space to template space via W=A'+B+C+D. Thus, only one
        interpolation of the epi data occurs.

    3) Will also censor volumes that have outlier, to not bias scaling
    """

    # Calculate T1-EPI rigid, T1-Template diffeo
    atlas_dir = "/home/data/madlab/atlases/vold2_mni"
    h_cmd = f"""
        cd {work_dir}

        align_epi_anat.py \
            -anat2epi \
            -anat struct+orig \
            -save_skullstrip \
            -suffix _al_junk \
            -epi epi_vrBase+orig \
            -epi_base 0 \
            -epi_strip 3dAutomask \
            -cost lpc+ZZ \
            -giant_move \
            -check_flip \
            -volreg off \
            -tshift off

        auto_warp.py -base {os.path.join(atlas_dir, "vold2_mni_brain+tlrc")} \
            -input struct_ns+orig -skull_strip_input no

        3dbucket -DAFNI_NIFTI_VIEW=tlrc -prefix struct_ns awpy/struct_ns.aw.nii*
        cp awpy/anat.un.aff.Xat.1D .
        cp awpy/anat.un.aff.qw_WARP.nii .
    """
    if not os.path.exists(os.path.join(work_dir, "struct_ns+tlrc.HEAD")):
        func_sbatch(h_cmd, 1, 4, 4, f"{subj_num}dif", work_dir)

    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        for run in epi_list:

            # determine input
            if blip_tog == 1:
                h_in = f"{run}_blip+orig"
            else:
                h_in = f"{run}+orig"

            # Calculate volreg for e/run
            h_cmd = f"""
                cd {work_dir}
                3dvolreg -verbose \
                    -zpad 1 \
                    -base epi_vrBase+orig \
                    -1Dfile dfile.{run}.1D \
                    -prefix {run}_volreg \
                    -cubic \
                    -1Dmatrix_save mat.{run}.vr.aff12.1D \
                    {h_in}
            """
            if not os.path.exists(os.path.join(work_dir, f"mat.{run}.vr.aff12.1D")):
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}vre", work_dir)

            # calc voxel dimension i
            h_cmd = (
                f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {run}+orig"
            )
            h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_gs_out = h_gs.communicate()[0]
            grid_size = h_gs_out.decode("utf-8").strip()

            # determine input
            if blip_tog == 1:
                h_nwarp = f"anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D blip_{run.split('_')[0]}_For_WARP+orig"
            else:
                h_nwarp = f"anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D"

            # concat warp calcs, apply to native EPI.
            #   Make a mask
            h_cmd = f"""
                cd {work_dir}

                cat_matvec -ONELINE \
                    anat.un.aff.Xat.1D \
                    struct_al_junk_mat.aff12.1D -I \
                    mat.{run}.vr.aff12.1D > mat.{run}.warp.aff12.1D

                3dNwarpApply -master struct_ns+tlrc \
                    -dxyz {grid_size} \
                    -source {run}+orig \
                    -nwarp '{h_nwarp}' \
                    -prefix {run}_warp

                3dcalc -overwrite -a {run}_blip+orig -expr 1 -prefix tmp_{run}_mask

                3dNwarpApply -master struct_ns+tlrc \
                    -dxyz {grid_size} \
                    -source tmp_{run}_mask+orig \
                    -nwarp 'anat.un.aff.qw_WARP.nii mat.{run}.warp.aff12.1D' \
                    -interp cubic \
                    -ainterp NN -quiet \
                    -prefix {run}_mask_warped

                3dTstat -min -prefix tmp_{run}_min {run}_mask_warped+tlrc
            """
            if not os.path.exists(os.path.join(work_dir, f"{run}_warp+tlrc.HEAD")):
                func_sbatch(h_cmd, 2, 4, 4, f"{subj_num}war", work_dir)

    # Determine minimum value, make mask
    #   beware the jabberwocky i.e. expanding braces in 3dMean
    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        h_cmd = f"""
            cd {work_dir}
            3dMean -datum short -prefix tmp_mean tmp_run-{{1..{len(epi_list)}}}_{phase}_min+tlrc
            3dcalc -a tmp_mean+tlrc -expr 'step(a-0.999)' -prefix {phase}_minVal_mask
        """
        if not os.path.exists(os.path.join(work_dir, f"{phase}_minVal_mask+tlrc.HEAD")):
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}min", work_dir)

        # make clean data
        epi_mask = os.path.join(work_dir, f"{phase}_minVal_mask+tlrc")
        for run in epi_list:
            h_warp = os.path.join(work_dir, f"{run}_warp+tlrc")
            h_clean = os.path.join(work_dir, f"{run}_volreg_clean")
            h_cmd = f"3dcalc -a {h_warp} -b {epi_mask} -expr 'a*b' -prefix {h_clean}"
            if not os.path.exists(h_clean + "+tlrc.HEAD"):
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}cle", work_dir)

    # %%
    """
    Step 5: Blur, Make masks

    1) Blur epi data - 1.5 * voxel dim, rounded up to nearest int. FWHM.

    2) Make a union mask, where sufficient signal exists for both T1w
        and T2*w at each voxel for analyses. Incorporated at the
        group-level analysis.

    3) Make tissue masks. The WM mask will be used later to derive
        nuissance regressors for the REML.
        Note: this references some custom masks, and is based in
            atropos rather than in AFNIs tiss seg protocol.
    """

    # Blur
    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        for run in epi_list:

            # calc voxel dim i
            h_cmd = (
                f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {run}+orig"
            )
            h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_gs_out = h_gs.communicate()[0]
            grid_size = h_gs_out.decode("utf-8").strip()
            blur_size = math.ceil(1.5 * float(grid_size))

            # do blur
            h_bin = os.path.join(work_dir, f"{run}_volreg_clean+tlrc")
            h_bout = os.path.join(work_dir, f"{run}_blur")
            h_cmd = f"3dmerge -1blur_fwhm {blur_size} -doall -prefix {h_bout} {h_bin}"
            if not os.path.exists(f"{h_bout}+tlrc.HEAD"):
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}blur", work_dir)

    # Make EPI-T1 union mask (mask_epi_anat)
    run_list = [
        x.split(".")[0]
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, "*blur+tlrc.HEAD")
    ]

    if not os.path.exists(os.path.join(work_dir, "mask_epi_anat+tlrc.HEAD")):

        for i in run_list:
            h_mout = os.path.join(work_dir, f"tmp_mask.{i}")
            h_min = os.path.join(work_dir, i)
            if not os.path.exists(h_mout + ".HEAD"):
                h_cmd = f"3dAutomask -prefix {h_mout} {h_min}"
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}mau", work_dir)

        h_cmd = f"""
            cd {work_dir}
            3dmask_tool -inputs tmp_mask.*+tlrc.HEAD -union -prefix tmp_mask_allRuns
            3dresample -master tmp_mask_allRuns+tlrc -input struct_ns+tlrc \
                -prefix tmp_anat_resamp
            3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_anat_resamp+tlrc \
                -prefix tmp_mask_struct
            3dmask_tool -input tmp_mask_allRuns+tlrc tmp_mask_struct+tlrc -inter \
                -prefix mask_epi_anat
            3dABoverlap -no_automask tmp_mask_allRuns+tlrc tmp_mask_struct+tlrc | \
                tee out.mask_ae_overlap.txt
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}uni", work_dir)

    # Make tissue-class masks
    #   I like Atropos better than AFNI's way, so use those priors
    atropos_dict = {1: "CSF", 2: "GMc", 3: "WM", 4: "GMs"}
    atropos_dir = os.path.join(atlas_dir, "priors_ACT")
    h_tcin = os.path.join(work_dir, run_list[0])

    for i in atropos_dict:
        h_tiss = atropos_dict[i]
        if not os.path.exists(
            os.path.join(work_dir, f"final_mask_{h_tiss}_eroded+tlrc.HEAD")
        ):
            h_cmd = f"""
                module load c3d/1.0.0
                cd {work_dir}

                c3d {atropos_dir}/Prior{i}.nii.gz -thresh 0.3 1 1 0 \
                    -o tmp_{h_tiss}_bin.nii.gz
                3dresample -master {h_tcin} -rmode NN -input tmp_{h_tiss}_bin.nii.gz \
                    -prefix final_mask_{h_tiss}+tlrc
                3dmask_tool -input tmp_{h_tiss}_bin.nii.gz -dilate_input -1 \
                    -prefix tmp_mask_{h_tiss}_eroded
                3dresample -master {h_tcin} -rmode NN -input tmp_mask_{h_tiss}_eroded+orig \
                    -prefix final_mask_{h_tiss}_eroded
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}atr", work_dir)

    # %%
    """
    Step 6: Scale data

    1) Data is scaled by mean signal
    """

    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        for run in epi_list:
            if not os.path.exists(os.path.join(work_dir, f"{run}_scale+tlrc.HEAD")):
                h_cmd = f"""
                    cd {work_dir}
                    3dTstat -prefix tmp_tstat_{run} {run}_blur+tlrc
                    3dcalc -a {run}_blur+tlrc \
                        -b tmp_tstat_{run}+tlrc \
                        -c {phase}_minVal_mask+tlrc \
                        -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
                        -prefix {run}_scale
                """
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}scale", work_dir)


# %%
def main():

    args = func_argparser().parse_args()
    h_data_dir = os.path.join(args.h_par, "dset", args.h_sub, args.h_ses)
    h_work_dir = os.path.join(args.h_par, "derivatives", args.h_sub, args.h_ses)

    if not os.path.exists(h_work_dir):
        os.makedirs(h_work_dir)

    # print(h_data_dir, h_work_dir, args.h_sub, args.h_ses, args.h_phl, args.h_blt)
    func_preproc(h_data_dir, h_work_dir, args.h_sub, args.h_ses, args.h_phl, args.h_blt)


if __name__ == "__main__":
    main()

# %%
