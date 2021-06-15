from os import system, getcwd, listdir, makedirs, rename, chdir
from os.path import splitext, isfile
from shutil import copyfile
import glob
from make_mdp import make_mdp
import numpy as np 
import matplotlib.pyplot as plt
from utils import plot_xvg, send_mail, gpu_manager
import pandas as pd
from multiprocessing import Pool
import argparse
from moles import init

init()
parser = argparse.ArgumentParser()
parser.add_argument("-gmx", type=str, help="gromacs --- default gmx", default = 'gmx')
parser.add_argument("-mdrun", type=str, help="mdrun --- default mdrun_tmpi", default = 'mdrun_tmpi')
parser.add_argument("-nt", "--numthread", type=int, help="number of threads --- default 4", default = 4)
parser.add_argument("-ns", "--nanoseconds", type=int, help="ns of simulation --- default 5 ns", default = 5)
parser.add_argument("-step", type=float, help="step in ps --- default 0.002", default = 0.002)
parser.add_argument("-vsite", type=str, choices=['hydrogens', 'aromatic'], help="vsite --- default None", default = '')
parser.add_argument("-native", action='store_true', help="Native Contact Analysis")
parser.add_argument("-p", "--pool", type=int, help="number of process --- default 1", default = 1)
parser.add_argument("-s", "--system", type=str, help="system PDB file")
parser.add_argument("-r", "--receptor", type=str, help="receptor pdb --- default protein.pdb", default = 'protein.pdb')
parser.add_argument("-f", "--report", type=str, help="path of report --- default $PWD/report.csv", default = '$PWD/report.csv')
parser.add_argument("-l", "--ligname", type=str, help="ligand common name or mol2 ligand file w/o ext --- default ZINC", default = 'ZINC')
args = parser.parse_args()

sysList = glob.glob(f'{args.ligname}*.mol2')

def main(system):

    filename, ext = splitext(system)
    makedirs(filename,  exist_ok=True)
    rename(system, f'{filename}/{mol2}')
    chdir(filename)

    if args.vsite: 
        system(f'{args.gmx} pdb2gmx -ff charmm36m -f {system} -vsite {args.vsite} -o {args.system}_gmx.pdb -water tip3p -ignh -p topol.top')
        print('IF YOU ARE USING VSITE, RISE THE STEP! 0.002 --> 0.004')
    else:           system(f'{args.gmx} pdb2gmx -ff charmm36m -f {args.system} -o {args.system}_gmx.pdb -water tip3p -ignh -p topol.top')

    system(f'{args.gmx} editconf -f {args.system}_gmx.pdb -o {system}_gmx.pdb -d 1.0 -quiet')
    system(f'{args.gmx} solvate -cp {args.system}_gmx.pdb -o {system}_gmx.pdb -p topol.top -quiet')

    ions_mdp = make_mdp('ions')
    system(f'{args.gmx} grompp -f {ions_mdp} -c {args.system}_gmx.pdb -p topol.top -o Complex_b4ion.tpr -maxwarn 10 -quiet')
    system(f'echo \'SOL\' | {args.gmx} -quiet genion -s Complex_b4ion.tpr -o Complex_4mini.pdb -neutral -conc 0.15 -p topol.top -quiet')

    #------------ MINIMIZATION
    deviceID = gpu_manager()
    mini_mdp = make_mdp(mdp = 'mini')
    system(f'{args.gmx} grompp -f {mini_mdp} -c Complex_4mini.pdb -r Complex_4mini.pdb -p topol.top -o mini.tpr -maxwarn 10')
    system(f'{args.mdrun} -deffnm mini -nt {args.numthread} -gpu_id {deviceID} -v')

    #------------- EQUILIBRATION
    equi_mdp = make_mdp(mdp = 'equi')
    system(f'{args.gmx} grompp -f {equi_mdp} -c mini.gro -r mini.gro -p topol.top -o equi.tpr -maxwarn 10')
    deviceID = gpu_manager()
    system(f'{args.mdrun} -deffnm equi -nt {args.numthread} -gpu_id {deviceID} -v')

    #------------- MD
    MD_mdp = make_mdp(mdp = 'MD', ns = args.nanoseconds, dt = args.step)
    system(f'{args.gmx} grompp -f {MD_mdp} -c equi.gro -p topol.top -o MD.tpr -maxwarn 10')
    deviceID = gpu_manager()
    system(f'{args.mdrun} -v -deffnm MD -nt {args.numthread} -gpu_id {deviceID}')
    chdir('../')

if __name__=='__main__':
    p = Pool(int(args.pool))
    p.map_async(main, ligList, chunksize=1).get()
    p.close()

