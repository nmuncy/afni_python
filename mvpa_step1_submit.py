# %%
import os
import fnmatch
import subprocess
from datetime import datetime
import time

# set up
deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
code_dir = "/home/nmuncy/compute/afni_python"
task_dict = {"loc": ["face", "scene", "num"]}
sess = "ses-S1"

# %%
"""
Step 1: keys
"""
subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
subj_list.sort()

# scan key
ref_subj = os.path.join(deriv_dir, subj_list[0])
h_cmd = f"module load afni-20.2.06 \n 3dinfo -tr {os.path.join(ref_subj, sess)}/run-1_{list(task_dict.keys())[0]}_scale+tlrc"
h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
len_tr = float(h_tr.communicate()[0].decode("utf-8").strip())
with open(os.path.join(deriv_dir, "scan_key.txt"), "w") as scan_key:
    scan_key.write(f"TR {len_tr}")

# study key
with open(os.path.join(deriv_dir, "study_key.txt"), "w") as study_key:
    study_key.write(f"{list(task_dict.keys())[0]} task")

# task, condition keys
if not os.path.exists(os.path.join(deriv_dir, "task_key.txt")):
    for count, task in enumerate(task_dict.keys()):

        # task key
        with open(os.path.join(deriv_dir, "task_key.txt"), "a") as task_key:
            task_key.write(f"task00{count+1} {task}")

        # condition keys
        cond_list = task_dict[task]
        for cc, cond in enumerate(cond_list):
            key_dir = os.path.join(deriv_dir, "models", f"model00{count+1}")
            if not os.path.exists(key_dir):
                os.makedirs(key_dir)
            with open(os.path.join(key_dir, "condition_key.txt"), "a") as cond_key:
                cond_key.write(f"task00{count+1} cond00{cc+1} {cond}\n")


# %%
"""
Step 2: submit
"""
# set up stdout/err capture
current_time = datetime.now()
out_dir = os.path.join(
    deriv_dir, f'Slurm_out/MVPA1_{current_time.strftime("%H%M_%d-%m-%y")}'
)
os.makedirs(out_dir)

for subj in subj_list:

    # Set stdout/err file
    h_out = os.path.join(out_dir, f"out_{subj}.txt")
    h_err = os.path.join(out_dir, f"err_{subj}.txt")

    # submit command
    sbatch_job = f"""
        sbatch \
        -J "PPI{i.split("-")[1]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
        -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
        python {code_dir}/mvpa_step1_setup.py {subj} {len_tr}"
    """

    sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_submit.communicate()[0]
    print(job_id)
    time.sleep(1)
