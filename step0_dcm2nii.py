# --- Notes:
#
# TODO: resolve the multiple T1w issue


import os
import subprocess
import tarfile
import fnmatch
import re
import time
import sys


# Submit jobs to slurm
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_sub, h_ses, h_str):

    full_name = h_sub + "_" + h_ses + "_" + h_str
    sbatch_job = f"sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --wrap='{command}'"
    # --account iacc_madlab --qos pq_madlab \

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
    sess = "ses-" + dcm_tar.split("-")[4].split(".")[0]
    subj = "sub-" + dcm_tar.split("-")[3]
    subj_dir = os.path.join(work_dir, subj, sess)
    dcm_dir = os.path.join(source_dir, subj, sess)

    for i in [subj_dir, dcm_dir]:
        if not os.path.exists(i):
            os.makedirs(i)

    tar_ball = os.path.join(data_dir, dcm_tar)
    tar_out = dcm_tar.split(".")[0]
    if not os.path.exists(os.path.join(dcm_dir, tar_out)):
        h_cmd = f"tar -C {dcm_dir} -xzf {tar_ball}"
        func_sbatch(h_cmd, 1, 1, 1, subj, sess, "tarEx")

    scan_list = os.listdir(os.path.join(dcm_dir, tar_out, "scans"))
    for i in scan_dict:

        h_list = [
            x for x in scan_list if i in x and "setter" not in x and "dMRI" not in x
        ]
        h_out_dir = os.path.join(subj_dir, scan_dict[i])
        if not os.path.exists(h_out_dir):
            os.makedirs(h_out_dir)

        for j in h_list:

            h_input_dir = os.path.join(
                dcm_dir, tar_out, "scans", j, "resources/DICOM/files"
            )

            # This might be cleaner in a class
            if scan_dict[i] == "func":
                h_out_str = (
                    subj
                    + "_"
                    + sess
                    + "_"
                    + "task-vCAT"
                    + "_run-"
                    + j.split("_")[-1]
                    + "_bold"
                )
            elif scan_dict[i] == "fmap":
                h_out_str = (
                    subj + "_" + sess + "_acq-func_dir-" + j.split("_")[-1] + "_epi"
                )
            elif scan_dict[i] == "anat":
                h_out_str = subj + "_" + sess + "_T1w"

            if not os.path.exists(os.path.join(h_out_dir, h_out_str)):
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
