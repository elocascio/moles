from os import system, getcwd, listdir, makedirs, rename, chdir
from os.path import splitext, isfile
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
import pandas as pd
import pickle as pkl
from rdkit.Chem import PandasTools, MolToSmiles, MolFromMol2File


name, lig_pattern, receptor, num = argv

ligList = glob.glob(f'{lig_pattern}*.mol2')

loop = 1
status_list = []
mean = []
smiles = []


print(receptor, len(ligList))

system(f'{gmx} pdb2gmx -ff charmm36m -f {receptor} -o Receptor_gmx.pdb -water tip3p -ignh -p topol.top {null}')

for mol2 in ligList:
    filename, ext = splitext(mol2)
    
    try:
        smile = MolToSmiles(MolFromMol2File(mol2))
    except:
        smile = 'C'
    ligpdb = filename + '.pdb'
    makedirs(filename,  exist_ok=True)
    rename(mol2, f'{filename}/{mol2}')
    chdir(filename)
    system(f'{MATCH} {mol2}')
    system(f'python {charmm2gmx} {getcwd()}/{filename}.rtf {getcwd()}/{filename}.prm {filename}.ff')
    ligand_ff = f'{filename}.ff'
    system(f'obabel -imol2 {mol2} -opdb -O {filename}.pdb')
    system(f'{gmx} pdb2gmx --ff {filename} -f {filename}.pdb -o Ligand_gmx.pdb -p Ligand.top {null}')

# WRITE ITP FILE
    if not isfile(f'{ligand_ff}/ffbonded.itp'): 
        status_list.append('ERROR')
        mean.append('ERROR')
        print('ERROR')
        sleep(1)
        continue
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
    
    system(f'{gmx} editconf -f {Complex.name} -o {Complex.name} -d 1.0 {null}')
    system(f'{gmx} solvate -cp {Complex.name} -o {Complex.name} -p topol.top {null}')

    ions_mdp = make_mdp('ions')
    system(f'{gmx} grompp -f {ions_mdp} -c {Complex.name} -p topol.top -o Complex_b4ion.tpr -maxwarn 10 {null}')
    system(f'echo 15 | {gmx} genion -s Complex_b4ion.tpr -o Complex_4mini.pdb -neutral -conc 0.15 -p topol.top {null}')

    system(f"""{gmx} make_ndx -f Complex_4mini.pdb -o index.ndx << EOF
1 | 20
q
EOF""")
    
#    system(f'echo "Protein_UNK" | {gmx} genrestr -f Complex_4mini.pdb -n index.ndx')
#    if os.path.isfile('')
    mini_mdp = make_mdp(mdp = 'mini')
    system(f'{gmx} grompp -f {mini_mdp} -c Complex_4mini.pdb -r Complex_4mini.pdb -p topol.top -o mini.tpr -maxwarn 10 {null}')
    system(f'{mdrun} -v -deffnm mini')
    
    equi_mdp = make_mdp(mdp = 'equi')
    system(f'{gmx} grompp -f {equi_mdp} -c mini.gro -r mini.gro -p topol.top -o equi.tpr -maxwarn 10 {null}')
    system(f'{mdrun} -v -deffnm equi')

    MD_mdp = make_mdp(mdp = 'MD')
    system(f'{gmx} grompp -f {MD_mdp} -c equi.gro -p topol.top -o MD.tpr -maxwarn 10 {null}')
    system(f'{mdrun} -v -deffnm MD')

    system(f'{gmx} editconf -f MD.tpr -o MD.pdb {null}')
    system(f"""{gmx} mindist -f MD.xtc -s MD.tpr -d 0.45 -n index.ndx -on contacts.xvg << EOF
    Protein
    13
    EOF""")

    time, contacts = plot_xvg('contacts.xvg', 'number' ,'time', 'concats')
    status, contacts_mean = detachmet(contacts)
    status_list.append(status)
    mean.append(contacts_mean)
    smiles.append(smile)

    pkl.dump(ligList, open('../ligList', 'wb'))
    pkl.dump(status_list, open('../status', 'wb'))
    pkl.dump(mean, open('../mean', 'wb'))
    
    if loop % int(num) == 0 or loop == len(ligList):
        df = pd.DataFrame({
            'ligand': ligList[:loop],
            'status': status_list,
            'contacts_average': mean,
            'smiles' : smiles})
        df = df[df['status'] != 'ERROR']
        df = df.sort_values(by=['contacts_average', 'status'], ascending = False)
        body = f"""
    This is an auto-generated email. Do not respond to this email address.
    
    This is {name} (aka {real_name}), reporting:

    HTMD
    Ligand processed : {loop}/{len(ligList)}
    Attached : {status_list.count('attached')}
    Detached : {status_list.count('detachment')}
    
    In allegato il report da aprire con browser!""" 
        PandasTools.AddMoleculeColumnToFrame(df, 'smiles', 'Molecule')
        report = df.to_html()
        with open('../Report.html', 'w') as html:
            html.write(report)

        send_mail(
            destination=['ettore.locascio@unicatt.it', 'alessandro.arcovito@unicatt.it'],
            subject = f"Report of Ligand-PALS1 {loop}/{len(ligList)}",
            content= body,
            attachment = "../Report.html",
            )
        
        df.to_pickle("../Report")
    else: pass
    loop += 1
    chdir('../')