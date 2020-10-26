"""
Notes:

Wrapper script for ppi_job.py.

Update paths in "set up" section.

decon_type can be dmBLOCK, GAM, or 2GAM
"""

# %%
import os
from datetime import datetime
import fnmatch
import subprocess
import json
import time


# set up
code_dir = "/home/nmuncy/compute/afni_python"
work_dir = "/scratch/madlab/nate_vCAT"
sess_list = ["ses-S1"]
phase_list = ["loc"]
decon_type = "2GAM"
seed_dict = {"LHC": "-24 -12 -22"}
stim_dur = 2


def main():

    # set up stdout/err capture
    deriv_dir = os.path.join(work_dir, "derivatives")
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, f'Slurm_out/PPI1_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    # submit job for each subj/sess/phase
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    for i in subj_list:
        for j in sess_list:
            subj_dir = os.path.join(deriv_dir, i, j)
            for k in phase_list:
                if not os.path.exists(
                    os.path.join(
                        subj_dir,
                        f"{k}_{decon_type}_ppi_stats_REML+tlrc.HEAD",
                    )
                ):

                    # write json to avoid quotation issue
                    with open(os.path.join(subj_dir, "seed_dict.json"), "w") as outfile:
                        json.dump(seed_dict, outfile)

                    # Set stdout/err file
                    h_out = os.path.join(out_dir, f"out_{i}_{j}_{k}.txt")
                    h_err = os.path.join(out_dir, f"err_{i}_{j}_{k}.txt")

                    # submit command
                    sbatch_job = f"""
                        sbatch \
                        -J "PPI{i.split("-")[1]}" -t 15:00:00 --mem=4000 --ntasks-per-node=1 \
                        -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
                        --account iacc_madlab --qos pq_madlab \
                        --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
                        python {code_dir}/ppi_step1_job.py {i} {j} {k} {decon_type} {deriv_dir} {stim_dur}"
                    """

                    sbatch_submit = subprocess.Popen(
                        sbatch_job, shell=True, stdout=subprocess.PIPE
                    )
                    job_id = sbatch_submit.communicate()[0]
                    print(job_id)
                    time.sleep(1)


if __name__ == "__main__":
    main()

# %%
