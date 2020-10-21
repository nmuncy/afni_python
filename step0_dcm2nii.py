"""
Notes:

1) Will extract tarball from data_dir to source_dir, then make
    NIfTI files in work_dir.

2) Will organize work_dir in BIDS format.

3) Gets submitted from wrapper script, will generate sbatch
    subprocesses for each dcm2nii command.

TODO:
  1) resolve the multiple T1w issue
  2) update to support multiple func phase names
        vCAT, loc are currently hardcoded output strings
"""
# %%
import os
import subprocess
import time
import sys


# Submit jobs to slurm, wait for job to finish
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_str, slurm_dir, subj, sess):

    full_name = f"{slurm_dir}/sbatch_writeOut_{subj}_{sess}_{h_str}"
    sbatch_job = f"""
        sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p centos7_IB_44C_512G  -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap="{command}"
    """
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_response.communicate()[0]
    print(job_id, h_str, sbatch_job)

    while_count = 0
    status = False
    while not status:

        check_cmd = "squeue -u $(whoami)"
        sq_check = subprocess.Popen(check_cmd, shell=True, stdout=subprocess.PIPE)
        out_lines = sq_check.communicate()[0]
        b_decode = out_lines.decode("utf-8")

        if h_str not in b_decode:
            status = True

        if not status:
            while_count += 1
            print(f"Wait count for sbatch job {h_str}: ", while_count)
            time.sleep(3)
    print(f'Sbatch job "{h_str}" finished')


def func_dcm2nii(scan_list, scan_name, scan_type, out_dir, scan_dir, subj_sess):

    for tmp_scan in scan_list:

        # Determine BIDS out string
        if scan_type == "func":
            run = tmp_scan.split("_")[-1]
            if scan_name == "Study":
                out_str = f"{subj_sess}_task-vCAT_run-{run}_bold"
            elif scan_name == "loc":
                out_str = f"{subj_sess}_task-loc_run-{run}_bold"
        elif scan_type == "fmap":
            blip = tmp_scan.split("_")[-1]
            out_str = f"{subj_sess}_acq-func_dir-{blip}_epi"
        elif scan_type == "anat":
            out_str = f"{subj_sess}_T1w"

        # write, submit job
        if not os.path.exists(os.path.join(out_dir, f"{out_str}.nii.gz")):
            h_cmd = f"""
                module load dcm2niix
                dcm2niix -b y -ba y -z y -o {out_dir} -f {out_str} {os.path.join(scan_dir, tmp_scan, "resources/DICOM/files")}
            """
            func_sbatch(h_cmd, 1, 1, 1, "dcm2nii", slurm_dir, subj, sess)


# %%
# def func_job(dcm_tar, data_dir, work_dir, source_dir, scan_dict, slurm_dir):

# for testing
dcm_tar = "Mattfeld_REVL-000-vCAT-005-S1.tar.gz"
data_dir = "/home/data/madlab/Mattfeld_vCAT/sourcedata"
work_dir = "/scratch/madlab/nate_vCAT/dset"
source_dir = "/scratch/madlab/nate_vCAT/sourcedata"
scan_dict = {"func": ["Study", "loc"], "anat": "T1w", "fmap": "Dist"}
slurm_dir = "/home/nmuncy/compute/afni_python"

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
    func_sbatch(h_cmd, 1, 1, 1, "tarEx", slurm_dir, subj, sess)

# %%
# make scans found in scan_dict
scan_dir = os.path.join(dcm_dir, tar_out, "scans")
for i in scan_dict:

    # make output dir
    out_dir = os.path.join(subj_dir, i)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # submit dcm2nii job
    if isinstance(scan_dict[i], list):
        for j in scan_dict[i]:
            h_list = [
                x
                for x in os.listdir(os.path.join(dcm_dir, tar_out, "scans"))
                if j in x and "setter" not in x and "dMRI" not in x and "32" not in x
            ]
            # print(h_list, j, i, out_dir, scan_dir, f"{subj}_{sess}")
            func_dcm2nii(h_list, j, i, out_dir, scan_dir, f"{subj}_{sess}")
    else:
        h_list = [
            x
            for x in os.listdir(os.path.join(dcm_dir, tar_out, "scans"))
            if scan_dict[i] in x and "setter" not in x and "dMRI" not in x
        ]
        # print(h_list, scan_dict[i], i, out_dir, scan_dir, f"{subj}_{sess}")
        func_dcm2nii(h_list, scan_dict[i], i, out_dir, scan_dir, f"{subj}_{sess}")


# %%
# def main():

#     h_dcm_tar = str(sys.argv[1])
#     h_data_dir = str(sys.argv[2])
#     h_par_dir = str(sys.argv[3])
#     h_slurm = str(sys.argv[4])

#     h_work_dir = os.path.join(h_par_dir, "dset")
#     h_source_dir = os.path.join(h_par_dir, "sourcedata")
#     h_scan_dict = {"Study": "func", "T1w": "anat", "Dist": "fmap", "loc": "func"}

#     func_job(h_dcm_tar, h_data_dir, h_work_dir, h_source_dir, h_scan_dict, h_slurm)


# if __name__ == "__main__":
#     main()
