"""
Notes:

Wrapper script for step1_preproc.py.

Update paths in "set up" section.

phase_list = list of phases gathered within a single session.
    For example, if a study and then a test phase were both scanned
    during the same session, then phase_list = ["study", "test"]

    Note: step1_preproc.py currently only accepts one phase per
        session, but this will be updated at some point.

TODO:
    1) combine sess/phase_list to dictionary?
"""

# %%
import os
from datetime import datetime
import subprocess
import time

# set up
code_dir = "/home/nmuncy/compute/afni_python"
work_dir = "/scratch/madlab/nate_vCAT"
sess_list = ["ses-S1"]
phase_list = ["vCAT"]
blip_toggle = 0  # 1 = on, 0 = off

# set up stdout/err capture
current_time = datetime.now()
out_dir = os.path.join(
    work_dir, f'derivatives/Slurm_out/TS1_{current_time.strftime("%H%M_%d-%m-%y")}'
)
os.makedirs(out_dir)

# submit job for each subj/sess/phase
subj_list = os.listdir(os.path.join(work_dir, "dset"))

# %%
for i in subj_list:
    for j in sess_list:
        for k in phase_list:
            if not os.path.exists(
                os.path.join(
                    work_dir, "derivatives", i, j, f"run-1_{k}_scale+tlrc.HEAD"
                )
            ):

                h_out = os.path.join(out_dir, f"out_{i}_{j}_{k}.txt")
                h_err = os.path.join(out_dir, f"err_{i}_{j}_{k}.txt")

                sbatch_job = f"""
                    sbatch \
                    -J "TS1" -t 10:00:00 --mem=4000 --ntasks-per-node=1 \
                    -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
                    --account iacc_madlab --qos pq_madlab \
                    --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
                    python {code_dir}/step1_preproc.py {i} {j} {k} {work_dir} {blip_toggle}"
                """

                sbatch_submit = subprocess.Popen(
                    sbatch_job, shell=True, stdout=subprocess.PIPE
                )
                job_id, error = sbatch_submit.communicate()
                print(job_id)

                time.sleep(1)
# %%
