# %%
import os
import subprocess
import tarfile
import fnmatch
import re
import time
import sys


# %%
test_mode = True
if test_mode:
    dcm_tar = "Mattfeld_REVL-000-vCAT-005-S1.tar.gz"
else:
    dcm_tar = str(sys.argv[1])


data_dir = "/home/data/madlab/Mattfeld_vCAT/sourcedata"
work_dir = "/scratch/madlab/nate_vCAT/dset"
scan_dict = {"REVL": "func", "T1w": "anat", "Dist": "fmap"}


def func_sbatch_nii(out_dir, file_str, dcm_dir, h_str, subj, sess):

    sbatch_job = f"sbatch \
        -J {h_str} -t 00:30:00 --mem=4000 --ntasks-per-node=4 \
        -p centos7_IB_44C_512G  -o dcm2nii_{subj}_{sess}_{h_str}.out \
        -e dcm2nii_{subj}_{sess}_{h_str}{h_str}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap='module load dcm2niix \n \
        dcm2niix -b y -ba y -z y -o {out_dir} -f {file_str} {dcm_dir}'"

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
            time.sleep(2)
    print('Sbatch jobs "{h_str}" finished')


# %%

# def func_ex_tar()

sess = "ses-" + dcm_tar.split("-")[4].split(".")[0]
subj = "sub-" + dcm_tar.split("-")[3]
subj_dir = os.path.join(work_dir, subj, sess)

if not os.path.exists(subj_dir):
    os.makedirs(subj_dir)


with tarfile.open(os.path.join(data_dir, dcm_tar), mode="r") as tar_file:
    # print(tar_file)
    tar_list = [x for x in tar_file.getnames() if "/resources" not in x]
    del tar_list[0]

    # for j in tar_list:
    j = tar_list[0]
    if (
        not re.search("dMRI", j)
        and not re.search("32ch", j)
        and not re.search("setter", j)
        and not re.search("ROI", j)
    ):
        if re.search("Dist", j):

            h_dir = "dir-" + j.split("_")[-1]
            file_str = subj + "_" + sess + "_acq-func_" + h_dir + "_epi"
            # dcm_dir = os.path.join(data_dir, j, "resources/DICOM/")
            # print(os.listdir(dcm_dir))
            
            for key in scan_dict:
                if key in j.split("/")[-1]


# def main()