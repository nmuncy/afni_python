"""
Notes:

Wrapper script for step0_dcm2nii.py.

Update paths in "set paths" section.
"""

# %%
import os
from datetime import datetime
import fnmatch
import subprocess
import time

# set paths
code_dir = "/home/nmuncy/compute/afni_python"
tar_dir = "/home/data/madlab/Mattfeld_vCAT/sourcedata"
work_dir = "/scratch/madlab/nate_vCAT"

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

# %%
for i in tar_list:

    tar_file = i.split("/")[-1]
    tar_str = tar_file.split(".")[0]
    h_out = os.path.join(slurm_dir, f"out_{tar_str}.txt")
    h_err = os.path.join(slurm_dir, f"err_{tar_str}.txt")

    sbatch_job = f"""
        sbatch \
        -J "TS0" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
        -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
        python {code_dir}/step0_dcm2nii.py {tar_file} {tar_dir} {work_dir} {slurm_dir}"
    """

    sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id, error = sbatch_submit.communicate()
    print(job_id)

    # give jobs time to start so we don't hit the "pending lockout"
    time.sleep(30)
# %%
