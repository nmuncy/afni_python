"""
Notes:

1) Will extract tarball from data_dir to source_dir, then make
    NIfTI files in work_dir.

2) Will organize work_dir in BIDS format.

3) Gets submitted from wrapper script, will generate sbatch
    subprocesses for each dcm2nii command.

TODO:
  1) resolve the multiple T1w issue
  2) update h_scan_dict[func] to hold task name, reflect in code
  3) update sbatch stdout/err write location
"""

import os
import subprocess
import tarfile
import fnmatch
import re
import time
import sys


# Submit jobs to slurm, wait for job to finish
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_sub, h_ses, h_str):

    full_name = h_sub + "_" + h_ses + "_" + h_str
    sbatch_job = f"sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap='{command}'"

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


def func_job(dcm_tar, data_dir, work_dir, source_dir, scan_dict):

    # get paths, make output dirs
    sess = "ses-" + dcm_tar.split("-")[4].split(".")[0]
    subj = "sub-" + dcm_tar.split("-")[3]
    subj_dir = os.path.join(work_dir, subj, sess)
    dcm_dir = os.path.join(source_dir, subj, sess)

    for i in [subj_dir, dcm_dir]:
        if not os.path.exists(i):
            os.makedirs(i)

    # unzip tarball
    tar_ball = os.path.join(data_dir, dcm_tar)
    tar_out = dcm_tar.split(".")[0]
    if not os.path.exists(os.path.join(dcm_dir, tar_out)):
        h_cmd = f"tar -C {dcm_dir} -xzf {tar_ball}"
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "tarEx")

    # make scans found in scan_dict
    scan_list = os.listdir(os.path.join(dcm_dir, tar_out, "scans"))
    for i in scan_dict:

        # exclude certain dirs
        h_list = [
            x for x in scan_list if i in x and "setter" not in x and "dMRI" not in x
        ]

        # make output dir
        h_out_dir = os.path.join(subj_dir, scan_dict[i])
        if not os.path.exists(h_out_dir):
            os.makedirs(h_out_dir)

        for j in h_list:

            # path to dcms
            h_input_dir = os.path.join(
                dcm_dir, tar_out, "scans", j, "resources/DICOM/files"
            )

            # Determine BIDS out string
            h_str = j.split("_")[-1]
            if scan_dict[i] == "func":
                h_out_str = f"{subj}_{sess}_task-vCAT_run-{h_str}_bold"
            elif scan_dict[i] == "fmap":
                h_out_str = f"{subj}_{sess}_acq-func_dir-{h_str}_epi"
            elif scan_dict[i] == "anat":
                h_out_str = f"{subj}_{sess}_T1w"

            # write, submit job
            if not os.path.exists(os.path.join(h_out_dir, f"{h_out_str}.nii.gz")):
                h_cmd = f"""
                    module load dcm2niix
                    dcm2niix -b y -ba y -z y -o {h_out_dir} -f {h_out_str} {h_input_dir}
                """
                func_sbatch(h_cmd, 1, 1, 1, subj, sess, "dcm2nii")


def main():
    # h_dcm_tar = "Mattfeld_REVL-000-vCAT-005-S1.tar.gz"
    h_dcm_tar = str(sys.argv[1])
    h_data_dir = "/home/data/madlab/Mattfeld_vCAT/sourcedata"
    h_work_dir = "/scratch/madlab/nate_vCAT/dset"
    h_source_dir = "/scratch/madlab/nate_vCAT/sourcedata"
    h_scan_dict = {"Study": "func", "T1w": "anat", "Dist": "fmap"}

    func_job(h_dcm_tar, h_data_dir, h_work_dir, h_source_dir, h_scan_dict)


if __name__ == "__main__":
    main()
