from os import system, getcwd, listdir, makedirs, rename, chdir
from os.path import splitext
from sys import argv
import glob
from var import *
from shutil import copyfile
from make_mdp import make_mdp
from time import sleep
import numpy as np 
import matplotlib.pyplot as plt
from utils import plot_xvg, send_mail, detachmet
from subprocess import Popen, PIPE

name, lig_pattern, receptor = argv

ligList = glob.glob(f'{lig_pattern}*.mol2')
loop = 1
status_list = []
print(receptor, len(ligList))

system(f'{gmx} pdb2gmx -ff charmm36m -f {receptor} -o Receptor_gmx.pdb -water tip3p -ignh -p topol.top')

for mol2 in ligList:
    filename, ext = splitext(mol2)
    ligpdb = filename + '.pdb'
    makedirs(filename,  exist_ok=True)
    rename(mol2, f'{filename}/{mol2}')
    chdir(filename)
    system(f'{MATCH} {mol2}')
    system(f'python {charmm2gmx} {getcwd()}/{filename}.rtf {getcwd()}/{filename}.prm {filename}.ff')
    ligand_ff = f'{filename}.ff'
    system(f'babel -imol2 {mol2} -opdb {filename}.pdb')
    system(f'{gmx} pdb2gmx --ff {filename} -f {filename}.pdb -o Ligand_gmx.pdb -p Ligand.top')


# WRITE ITP FILE

    with open(f'{ligand_ff}/ffbonded.itp') as ffbonded, open(f'{ligand_ff}/ffnonbonded.itp') as ffnonbonded, open(f'Ligand.top') as ligand_top, open(f'{filename}.itp', 'w') as ligand_itp:
        ffnonbonded_lines = ffnonbonded.readlines()
        ffbonded_lines = ffbonded.readlines()
        ffnonbonded_lines += ffbonded_lines
        ligand_top_lines = ligand_top.readlines()
        ligand_itp.writelines(ffnonbonded_lines)
        for line in ligand_top_lines:
            if line.startswith(';') or line.startswith('#'):
                pass
            elif 'Other' in line:
                line = 'UNK   3'
                ligand_itp.write(line)
            elif '[ system ]' in line:
                break
            else: ligand_itp.write(line)
        print(f'{ligand_itp.name} scritto!')

# COPY AND MODIFY TOPOL.TOP
    copyfile('../topol.top', './topol.top')
    topol_line = open('topol.top').readlines()
    with open('topol.top', 'w') as topol:
        for line in topol_line:
            topol.write(line)
            if 'forcefield.itp' in line:
                topol.write(f'#include "{ligand_itp.name}"\n')
        topol.write('UNK   1\n')

# MAKE COMPLEX PDB

    with open('../Receptor_gmx.pdb') as receptor_gmx, open(f'{filename}.pdb') as ligand_pdb, open('Complex.pdb', 'w') as Complex:
        for line in receptor_gmx.readlines():
            if line.startswith('ATOM'):
                Complex.write(line)
        for line in ligand_pdb.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                Complex.write(line)
    
    system(f'{gmx} editconf -f {Complex.name} -o {Complex.name} -d 1.0')
    system(f'{gmx} solvate -cp {Complex.name} -o {Complex.name} -p topol.top')

    ions_mdp = make_mdp('ions')
    system(f'{gmx} grompp -f {ions_mdp} -c {Complex.name} -p topol.top -o Complex_b4ion.tpr -maxwarn 10')
    system(f'echo 15 | {gmx} genion -s Complex_b4ion.tpr -o Complex_4mini.pdb -neutral -conc 0.15 -p topol.top')

    system(f"""{gmx} make_ndx -f Complex_4mini.pdb -o index.ndx << EOF
1 | 20
q
EOF""")
    
#    system(f'echo "Protein_UNK" | {gmx} genrestr -f Complex_4mini.pdb -n index.ndx')
#    if os.path.isfile('')
    mini_mdp = make_mdp(mdp = 'mini')
    system(f'{gmx} grompp -f {mini_mdp} -c Complex_4mini.pdb -r Complex_4mini.pdb -p topol.top -o mini.tpr -maxwarn 10')
    system(f'{mdrun} -v -deffnm mini')
    
    equi_mdp = make_mdp(mdp = 'equi')
    system(f'{gmx} grompp -f {equi_mdp} -c mini.gro -r mini.gro -p topol.top -o equi.tpr -maxwarn 10')
    system(f'{mdrun} -v -deffnm equi')

    MD_mdp = make_mdp(mdp = 'MD')
    system(f'{gmx} grompp -f {MD_mdp} -c equi.gro -p topol.top -o MD.tpr -maxwarn 10')
    system(f'{mdrun} -v -deffnm MD')

    system(f'{gmx} editconf -f MD.tpr -o MD.pdb')
    system(f"""{gmx} mindist -f MD.xtc -s MD.tpr -d 0.45 -n index.ndx -on contacts.xvg << EOF
    Protein
    13
    EOF""")

    time, contacts = plot_xvg('contacts.xvg', 'number' ,'time', 'concats')
    status = detachmet(contacts)
    status_list.append(status)

    loop += 1
    if loop % 5 == 0:
        content = f"""
    This is an auto-generated email. Do not respond to this email address.
    
    HTMD
    Ligand processed : {loop}/{len(ligList)}
    Attached : {status_list.count('attached')}
    Detached : {status_list.count('detachment')}""" 

        send_mail('ettore.locascio@unicatt.it', content)
    else: pass

    chdir('../')