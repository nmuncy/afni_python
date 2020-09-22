# %%
"""
Notes

1) Only written to accept one experiment phase (e.g. Study, Test) per session, assumes
    all epi scans in session pertain to phase

2) Script will generate named sbatch subprocesses named for each sub-step.
    Output for each sbatch step is written to work_dir and prepeneded
    with "sbatch_writeOut"

3) Pre-processing steps include copying data to derivatives, distortion correcting
    volume registration, normalization, generating intersection masks, and
    scaling the data.

4) Written in Python 3.8, has afni and c3d dependencies.


TODO:
  1) Update template? Referencing special template/priors
  2) Update for multiple phases per experiment
  3) Receive work_dir from wrapper script
  4) Use something besides epi_dict?
        Originally I was going to pull info from the json sidecars, but they
        just don't contain the info I want.
"""

import json
import os
import sys
import subprocess
import fnmatch
import math
import time


# Submit jobs to slurm, check & wait for job to finish
#   Note: len(h_str) < 8
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_str, work_dir):

    full_name = f"{work_dir}/sbatch_writeOut_{h_str}"
    sbatch_job = f"""
        sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load afni-20.2.06 \n {command}"
    """
    print(h_str, sbatch_job)
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


# %%
# Do general pre-processing steps
def func_preproc(data_dir, work_dir, subj, sess, phase):

    """
    Step 1: Copy data into work_dir

    Get func, anat, fmap data. Rename appropriately.

    To account for different num of fmaps, will
    produce AP, PA fmap per run.

    Make epi_dict for later use, like pulling header info
    or looping through runs.
    """

    # # For testing
    # subj = "sub-005"
    # sess = "ses-S1"
    # phase = "vCAT"

    # par_dir = "/scratch/madlab/nate_vCAT"
    # data_dir = os.path.join(par_dir, "dset", subj, sess)
    # work_dir = os.path.join(par_dir, "derivatives", subj, sess)

    # if not os.path.exists(work_dir):
    #     os.makedirs(work_dir)

    # struct
    struct_nii = os.path.join(data_dir, "anat", "{}_{}_T1w.nii.gz".format(subj, sess))
    struct_raw = os.path.join(work_dir, "struct+orig")
    if not os.path.exists(struct_raw + ".HEAD"):
        h_cmd = "3dcopy {} {}".format(struct_nii, struct_raw)
        func_sbatch(h_cmd, 1, 1, 1, "struct", work_dir)

    # epi - only task, not rsFMRI
    epi_list = [
        epi
        for epi in os.listdir(os.path.join(data_dir, "func"))
        if fnmatch.fnmatch(epi, f"*task-{phase}*.nii.gz")
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
            func_sbatch(h_cmd, 1, 1, 1, "epi", work_dir)

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
            func_sbatch(h_cmd, 1, 1, 1, "fmap", work_dir)

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

    # step check
    if not os.path.exists(os.path.join(work_dir, f"run-1_{phase}+orig.HEAD")):
        print("Check scan missing. Exiting.")
        exit

    # %%
    """
    Step 2: Detect outliers voxels, blip correct

    1) Determine the proportion of voxels per volume that have outlier signal.
        Censor volumes that exceed limit.

    2) Correct for signal fallout using fmap. This approach is taken from
        afni_proc. It uses the fmap to "unwarp" the run epi.
    """

    for i in epi_dict.keys():

        # determine polort arg, file strings, find outliers
        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -tr {i}+orig"
        h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_len_tr = h_tr.communicate()[0]
        len_tr = h_len_tr.decode("utf-8").strip()

        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -ntimes {i}+orig"
        h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_nvol_out = h_nvol.communicate()[0]
        num_tr = h_nvol_out.decode("utf-8").strip()

        pol = 1 + math.ceil((float(len_tr) * float(num_tr)) / 150)

        if not os.path.exists(os.path.join(work_dir, f"out.cen.{i}.1D")):
            h_cmd = f"""
                cd {work_dir}
                3dToutcount -automask -fraction -polort {pol} \
                    -legendre {i}+orig > outcount.{i}.1D
                1deval -a outcount.{i}.1D -expr '1-step(a-0.1)' > out.cen.{i}.1D
            """
            func_sbatch(h_cmd, 1, 1, 1, "outlier", work_dir)

        h_run = i.split("_")[0]
        for j in ["Forward", "Reverse"]:

            # create median datasets and masks
            if not os.path.exists(
                os.path.join(work_dir, f"tmp_blip_med_masked_{h_run}_{j}+orig.HEAD")
            ):
                h_cmd = f"""
                    cd {work_dir}
                    3dTstat -median -prefix tmp_blip_med_{h_run}_{j} blip_{h_run}_{j}+orig
                    3dAutomask -apply_prefix tmp_blip_med_masked_{h_run}_{j} \
                        tmp_blip_med_{h_run}_{j}+orig
                """
                func_sbatch(h_cmd, 1, 1, 1, "fmap_med", work_dir)

        # comput midpoint warp, unwarp run data (solve for fall out), apply header
        if not os.path.exists(os.path.join(work_dir, f"{i}_blip+orig.HEAD")):
            h_cmd = f"""
                cd {work_dir}

                3dQwarp -plusminus -pmNAMES Rev For \
                    -pblur 0.05 0.05 -blur -1 -1 \
                    -noweight -minpatch 9 \
                    -source tmp_blip_med_masked_{h_run}_Reverse+orig \
                    -base tmp_blip_med_masked_{h_run}_Forward+orig \
                    -prefix blip_{h_run}

                3dNwarpApply -quintic -nwarp blip_{h_run}_For_WARP+orig \
                    -source {i}+orig -prefix {i}_blip

                3drefit -atrcopy blip_{h_run}_Forward+orig IJK_TO_DICOM_REAL {i}_blip+orig
            """
            func_sbatch(h_cmd, 1, 4, 2, "qwarp", work_dir)

    # step check
    for i in epi_dict.keys():
        if not os.path.exists(os.path.join(work_dir, f"{i}_blip+orig.HEAD")):
            print(f"File {i} missing. Exiting.")
            exit

    # %%
    """
    Step 3: Make volreg base

    1) The volreg base (epi_vr_base) is the single volume in entire experiment
    phase with the smallest number of outlier volumes.

    Note: If data is not time shifted, do so before determining volreg_base.
        Also, it was easier to write a small bash script due to the multiple
        afni commands.
    """

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

                3dbucket -prefix epi_vr_base ${{baseRun}}_blip+orig"[${{minouttr}}]"
                """.format(
                    work_dir
                )
            )

        h_cmd = f"source {vr_script}"
        func_sbatch(h_cmd, 1, 1, 1, "vrbase", work_dir)

    # step check
    if not os.path.exists(os.path.join(work_dir, "epi_vr_base+orig.HEAD")):
        print("VR base missing. Exiting.")
        exit

    # %%
    """
    Step 4: Calc, Perfrom normalization

    This step will perform the rigid alignments of T1-EPI (A)
    EPI-EPI base volume (B), and non-linear diffeomorphic of T1-Template (C).

    It will then concatenate these warp matrices, and warp EPI data from
    raw/native space to template space via W=A'+B+C. Thus, only one
    interpolation of the epi data occurs.

    Will also censor volumes that have outlier, to not bias scaling
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
            -epi epi_vr_base+orig \
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
        func_sbatch(h_cmd, 4, 4, 4, "diffeo", work_dir)

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
            func_sbatch(h_cmd, 1, 1, 1, "volreg", work_dir)

    # Concat calcs, warp EPI. Make mask.
    #   Could be combined with loop above
    for i in epi_dict.keys():

        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {i}+orig"
        h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_gs_out = h_gs.communicate()[0]
        grid_size = h_gs_out.decode("utf-8").strip()

        h_cmd = f"""
            cd {work_dir}

            cat_matvec -ONELINE \
                anat.un.aff.Xat.1D \
                struct_al_junk_mat.aff12.1D -I \
                mat.{i}.vr.aff12.1D > mat.{i}.warp.aff12.1D

            3dNwarpApply -master struct_ns+tlrc \
                -dxyz {grid_size} \
                -source {i}+orig \
                -nwarp 'anat.un.aff.qw_WARP.nii mat.{i}.warp.aff12.1D blip_{i}_For_WARP+orig' \
                -prefix {i}_warp

            3dcalc -overwrite -a {i}_blip+orig -expr 1 -prefix tmp_{i}_mask

            3dNwarpApply -master struct_ns+tlrc \
                -dxyz {grid_size} \
                -source tmp_{i}_mask+orig \
                -nwarp 'anat.un.aff.qw_WARP.nii mat.{i}.warp.aff12.1D' \
                -interp cubic \
                -ainterp NN -quiet \
                -prefix {i}_mask_warped

            3dTstat -min -prefix tmp_{i}_min {i}_mask_warped+tlrc
        """
        if not os.path.exists(os.path.join(work_dir, f"{i}_warp+tlrc.HEAD")):
            func_sbatch(h_cmd, 1, 4, 4, "warp", work_dir)

    # Determine minimum value, make mask
    #   wrote expanding braces into 3dmean command
    h_cmd = f"""
        cd {work_dir}
        3dMean -datum short -prefix tmp_mean tmp_run-{{1..{len(epi_dict.keys())}}}_{phase}_min+tlrc
        3dcalc -a tmp_mean+tlrc -expr 'step(a-0.999)' -prefix {phase}_minVal_mask
    """
    if not os.path.exists(os.path.join(work_dir, f"{phase}_minVal_mask+tlrc.HEAD")):
        func_sbatch(h_cmd, 1, 1, 1, "minVal", work_dir)

    # make clean data
    epi_mask = os.path.join(work_dir, f"{phase}_minVal_mask+tlrc")

    for i in epi_dict.keys():
        h_warp = os.path.join(work_dir, f"{i}_warp+tlrc")
        h_clean = os.path.join(work_dir, f"{i}_volreg_clean")
        h_cmd = f"3dcalc -a {h_warp} -b {epi_mask} -expr 'a*b' -prefix {h_clean}"

        if not os.path.exists(h_clean + "+tlrc.HEAD"):
            func_sbatch(h_cmd, 1, 1, 1, "clean", work_dir)

    # step check
    for i in epi_dict.keys():
        if not os.path.exists(os.path.join(work_dir, f"{i}_volreg_clean+tlrc.HEAD")):
            print(f"{i}_volreg_clean+tlrc missing. Exiting.")
            exit

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
    for i in epi_dict.keys():

        h_cmd = f"module load afni-20.2.06 \n cd {work_dir} \n 3dinfo -di {i}+orig"
        h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_gs_out = h_gs.communicate()[0]
        grid_size = h_gs_out.decode("utf-8").strip()
        blur_size = math.ceil(1.5 * float(grid_size))

        h_bin = os.path.join(work_dir, f"{i}_volreg_clean+tlrc")
        h_bout = os.path.join(work_dir, f"{i}_blur")
        h_cmd = f"3dmerge -1blur_fwhm {blur_size} -doall -prefix {h_bout} {h_bin}"
        if not os.path.exists(f"{h_bout}+tlrc.HEAD"):
            func_sbatch(h_cmd, 1, 1, 1, "blur", work_dir)

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
                func_sbatch(h_cmd, 1, 1, 1, "mauto", work_dir)

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
        func_sbatch(h_cmd, 1, 1, 1, "union", work_dir)

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
            func_sbatch(h_cmd, 1, 1, 1, "atropos", work_dir)

    # step check
    if not os.path.exists(os.path.join(work_dir, "final_mask_WM_eroded+tlrc.HEAD")):
        print("final_mask_WM_eroded+tlrc missing. Exiting.")
        exit

    # %%
    """
    Step 6: Scale data

    Data is scaled by mean signal
    """

    for i in epi_dict.keys():
        if not os.path.exists(os.path.join(work_dir, f"{i}_scale+tlrc.HEAD")):
            h_cmd = f"""
                cd {work_dir}
                3dTstat -prefix tmp_tstat_{i} {i}_blur+tlrc
                3dcalc -a {i}_blur+tlrc \
                    -b tmp_tstat_{i}+tlrc \
                    -c {epi_mask} \
                    -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
                    -prefix {i}_scale
            """
            func_sbatch(h_cmd, 1, 1, 1, "scale", work_dir)

    # step check
    for i in epi_dict.keys():
        if not os.path.exists(os.path.join(work_dir, f"{i}_scale+tlrc.HEAD")):
            print(f"{i}_scale+tlrc missing. Exiting.")
            exit


# %%
def main():

    subj = str(sys.argv[1])
    sess = str(sys.argv[2])
    phase = str(sys.argv[3])

    par_dir = "/scratch/madlab/nate_vCAT"
    data_dir = os.path.join(par_dir, "dset", subj, sess)
    work_dir = os.path.join(par_dir, "derivatives", subj, sess)

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    func_preproc(data_dir, work_dir, subj, sess, phase)


if __name__ == "__main__":
    main()

# %%
