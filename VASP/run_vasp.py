import logging                                                                               
logging.basicConfig(filename=f'run_vasp.log', filemode='a',
        level=logging.INFO, format="%(asctime)-15s %(message)s")          
logging.getLogger().addHandler(logging.StreamHandler())

import numpy as np
from os import mkdir, chdir
from os.path import isfile, isdir
from shutil import copy
from subprocess import call

command=os.environ.get("RUN_VASP", "srun vasp_std")
## debug
# command = "cp POSCAR CONTCAR"
output="vasp.out"

if isfile("previous_run"):
    previous_run = np.loadtxt("previous_run")
    previous_run = previous_run.reshape([-1])
    if len(previous_run) > 0:
        folder_name = int(previous_run[-1])
    else:
        previous_run = []
        folder_name = 0
else:
    previous_run = []
    folder_name = 0

logging.info(f"start folder_name {folder_name}")
while isdir(f"{folder_name}"):
    if isfile(f"{folder_name}/CONTCAR"):
        copy(f"{folder_name}/CONTCAR", "POSCAR")
        logging.info(f"cp {folder_name}/CONTCAR POSCAR")
    else:
        raise RuntimeError(f"previous folder {folder_name} does not have legit CONTCAR")
    folder_name += 1
mkdir(f"{folder_name}")
logging.info(f"mkdir {folder_name}")
previous_run = np.hstack([previous_run, [folder_name]])
np.savetxt("previous_run", previous_run)

for filename in ["INCAR", "POSCAR", "KPOINTS", "POTCAR"]:
    copy(filename, f"{folder_name}")
    logging.info(f"cp {filename} {folder_name}")

chdir(f"{folder_name}")
logging.info(f"cd {folder_name}")
with open(output, "w+") as fout:
    logging.info(f"{command} > {output}")
    call(command.split(), stdout=fout)
