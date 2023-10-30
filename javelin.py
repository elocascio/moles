import argparse
from os import getcwd
from os.path import isfile
import subprocess
from time import sleep
from os import system

parser = argparse.ArgumentParser()
parser.add_argument("-name", type=str, help="job name")
parser.add_argument("-nt", type=int, help="cores", default=24)
parser.add_argument("-user", type=str, help="username", default="elocasci")
parser.add_argument("-email", type=str, help="email", default="elocascio@unicatt.it")
args = parser.parse_args()

folder = getcwd()
print(folder)
cpt = ""


while True:
    qme = (
        (
            subprocess.Popen(
                ["squeue", "-u", args.user], stdout=subprocess.PIPE
            ).communicate()[0]
        ).decode("ascii")
    ).split()
    if args.name in qme:
        print(f"found {args.name} in {qme}")
        sleep(120)
        continue
    else:
        with open("m100.sh", "w") as sh:
            if isfile("MD.cpt"):
                cpt = "-cpi"
                print("checkpoint found!")

            sh.write(
                f"""#!/bin/bash
#SBATCH -A IscrC_INSIDE1
#SBATCH -p m100_usr_prod
#SBATCH --qos=m100_qos_dbg
#SBATCH --time 2:00:00     # format: HH:MM:SS
#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node={args.nt} # 8 tasks out of 128
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1        # 1 gpus per node out of 4
#SBATCH --mem=7100          # memory per node out of 246000MB
#SBATCH --job-name={args.name}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=args.email

module load profile/lifesc

module load autoload gromacs/2021.4

cd {folder}
export  OMP_NUM_THREADS={args.nt}
mdrun_thread_mpi -v -deffnm MD {cpt}
"""
            )
    system("sbatch m100.sh")
    sleep(1800)
