"""
Notes:

Wrapper script for step3_decon.py.

Update paths in "set up" section.

decon_type can be dmBLOCK, GAM, 2GAM, or TENT
"""

# %%
import os
from datetime import datetime
import fnmatch
import subprocess

# set up
code_dir = "/home/nmuncy/compute/afni_python"
work_dir = "/scratch/madlab/nate_vCAT"
sess_list = ["ses-S1"]
phase_list = ["vCAT"]
decon_type = "TENT"

# set up stdout/err capture
deriv_dir = os.path.join(work_dir, "derivatives")
current_time = datetime.now()
out_dir = os.path.join(
    deriv_dir, f'Slurm_out/TS3_{current_time.strftime("%H%M_%d-%m-%y")}'
)
os.makedirs(out_dir)

# submit job for each subj/sess/phase
subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]

# %%
for i in subj_list:
    for j in sess_list:
        for k in phase_list:
            if not os.path.exists(
                os.path.join(
                    deriv_dir,
                    i,
                    j,
                    f"{k}_{decon_type}_stats_REML+tlrc.HEAD",
                )
            ):

                h_out = os.path.join(out_dir, f"out_{i}_{j}_{k}.txt")
                h_err = os.path.join(out_dir, f"err_{i}_{j}_{k}.txt")

                sbatch_job = f"""
                    sbatch \
                    -J "TS3" -t 10:00:00 --mem=4000 --ntasks-per-node=1 \
                    -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
                    --account iacc_madlab --qos pq_madlab \
                    --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
                    python {code_dir}/step3_decon.py {i} {j} {k} {decon_type} {deriv_dir}"
                """

                sbatch_submit = subprocess.Popen(
                    sbatch_job, shell=True, stdout=subprocess.PIPE
                )
                job_id, error = sbatch_submit.communicate()
                print(job_id)
# %%
