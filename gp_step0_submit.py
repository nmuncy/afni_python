"""
Notes:

Wrapper script for step0_dcm2nii.py.

Usage - update "set up" section.
"""

# %%
import os
from datetime import datetime
import fnmatch
import subprocess
import time
import json

# set up
code_dir = "/home/nmuncy/compute/afni_python"
tar_dir = "/home/data/madlab/Mattfeld_vCAT/sourcedata"
work_dir = "/scratch/madlab/nate_vCAT"
scan_dict = {"func": ["Study", "loc"], "anat": "T1w", "fmap": "Dist"}


# %%
def main():

    # set up out_dir to capture stdout/err
    current_time = datetime.now()
    out_dir = f'derivatives/Slurm_out/TS0_{current_time.strftime("%H%M_%d-%m-%y")}'
    slurm_dir = os.path.join(work_dir, out_dir)

    # set up work_dir
    dir_list = ["dset", "sourcedata", "derivatives", out_dir]
    for i in dir_list:
        h_dir = os.path.join(work_dir, i)
        if not os.path.exists(h_dir):
            os.makedirs(h_dir)

    # list of tar balls
    tar_list = [x for x in os.listdir(tar_dir) if fnmatch.fnmatch(x, "*.tar.gz")]
    tar_list.sort()

    # write json to avoid quotation issue
    with open(os.path.join(slurm_dir, "scan_dict.json"), "w") as outfile:
        json.dump(scan_dict, outfile)

    # submit jobs
    for i in tar_list:
        # i = tar_list[0]

        tar_file = i.split("/")[-1]
        tar_str = tar_file.split(".")[0]
        h_out = os.path.join(slurm_dir, f"out_{tar_str}.txt")
        h_err = os.path.join(slurm_dir, f"err_{tar_str}.txt")

        sbatch_job = f"""
            sbatch \
            -J "GP0{i.split("-")[3]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
            -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
            python {code_dir}/gp_step0_dcm2nii.py {tar_file} {tar_dir} {work_dir} {slurm_dir}"
        """

        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)

        time.sleep(1)


if __name__ == "__main__":
    main()
# %%
