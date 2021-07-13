from os import system, getcwd, makedirs, rename, chdir
from os.path import splitext, isfile
import glob
from var import *
from shutil import copyfile
from make_mdp import make_mdp
import numpy as np 
from utils import plot_xvg, send_mail, detachmet, gpu_manager
import pandas as pd
from rdkit.Chem import MolToSmiles, MolFromMol2File, MolFromSmiles
from multiprocessing import Pool
import platform
from plip_contacts import *
import argparse
from Misc.moles import init

init()

parser = argparse.ArgumentParser()
parser.add_argument("-gmx", type=str, help="gromacs --- default gmx", default = 'gmx')
parser.add_argument("-mdrun", type=str, help="mdrun --- default mdrun_tmpi", default = 'mdrun_tmpi')
parser.add_argument("-nt", "--numthread", type=int, help="number of threads --- default 4", default = 4)
parser.add_argument("-ns", "--nanoseconds", type=int, help="ns of simulation --- default 5 ns", default = 5)
parser.add_argument("-step", type=float, help="step in ps --- default 0.002", default = 0.002)
parser.add_argument("-ln", "--lignum", type=int, help="send mail every n ligands --- default 50", default = 50)
parser.add_argument("-p", "--pool", type=int, help="number of process --- default 4", default = 4)
parser.add_argument("-l", "--ligname", type=str, help="ligand common name or mol2 ligand file w/o ext --- default ZINC", default = 'ZINC')
parser.add_argument("-i", "--ligID", type=str, help="ligand id --- default UNK", default = 'UNK')
parser.add_argument("-r", "--receptor", type=str, help="receptor pdb --- default protein.pdb", default = 'protein.pdb')
parser.add_argument("-f", "--report", type=str, help="path of report --- default ./report.csv", default = './report.csv')
args = parser.parse_args()

ligList = glob.glob(f'{args.ligname}*.mol2')

print(args.receptor, len(ligList))

system(f'{args.gmx} pdb2gmx -ff charmm36m -f {args.receptor} -o Receptor_gmx.pdb -water tip3p -ignh -p topol.top {null}')

def main(mol2):

    print(getcwd())
    if not isfile(mol2):
        return print('ERROR no mol2')
    
    filename, ext = splitext(mol2)
    
    try:
        smile = MolToSmiles(MolFromMol2File(mol2))
    except:
        smile = 'C'

    makedirs(filename,  exist_ok=True)
    rename(mol2, f'{filename}/{mol2}')
    chdir(filename)
    system(f'{MATCH} {mol2} {null}')
    if not isfile(f'{filename}.prm'):
        chdir('../')
        return print('ERROR')
        
    system(f'python {charmm2gmx} {getcwd()}/{filename}.rtf {getcwd()}/{filename}.prm {filename}.ff')
    ligand_ff = f'{filename}.ff'
    system(f'obabel -imol2 {mol2} -opdb -O {filename}.pdb {null}')
    system(f'{args.gmx} pdb2gmx --ff {filename} -f {filename}.pdb -o Ligand_gmx.pdb -p Ligand.top {null}')

# WRITE ITP FILE
    if not isfile(f'{ligand_ff}/ffbonded.itp'):
        chdir('../')
        return print('ERROR')
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
        ligand_itp.write("""
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif""")
        print(f'{ligand_itp.name} scritto!')

# COPY AND MODIFY TOPOL.TOP
    copyfile('../topol.top', './topol.top')
    copyfile('../posre.itp', './protein.itp')
    topol_line = open('topol.top').readlines()

    with open('topol.top', 'w') as topol:
        for line in topol_line:
            if 'forcefield.itp' in line:
                topol.write(line)
                topol.write(f'#include "{ligand_itp.name}"\n')
            elif 'posre.itp' in line:
                line = '#include "protein.itp"\n'
                topol.write(line)
            else: topol.write(line)
        topol.write('UNK   1\n')

# MAKE COMPLEX PDB

    with open('../Receptor_gmx.pdb') as receptor_gmx, open(f'{filename}.pdb') as ligand_pdb, open('Complex.pdb', 'w') as Complex:
        for line in receptor_gmx.readlines():
            if line.startswith('ATOM'):
                Complex.write(line)
        for line in ligand_pdb.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                Complex.write(line)
    
    system(f'{args.gmx} editconf -f {Complex.name} -o {Complex.name} -d 1.0 -quiet')
    system(f'{args.gmx} solvate -cp {Complex.name} -o {Complex.name} -p topol.top -quiet')

    ions_mdp = make_mdp('ions')
    system(f'{args.gmx} grompp -f {ions_mdp} -c {Complex.name} -p topol.top -o Complex_b4ion.tpr -maxwarn 10 -quiet {null}')
    system(f'echo 15 | {args.gmx} -quiet genion -s Complex_b4ion.tpr -o Complex_4mini.pdb -neutral -conc 0.15 -p topol.top -quiet {null}')

    system(f"""{args.gmx} -quiet make_ndx -f Complex_4mini.pdb -o index.ndx << EOF
1 | 20
q
EOF""")

#------------ MINIMIZATION
    deviceID = gpu_manager()
    mini_mdp = make_mdp(mdp = 'mini')
    system(f'{args.gmx} grompp -f {mini_mdp} -c Complex_4mini.pdb -r Complex_4mini.pdb -p topol.top -o mini.tpr -maxwarn 10')
    system(f'{args.mdrun} -deffnm mini -nt {args.numthread} -gpu_id {deviceID}')

#------------- EQUILIBRATION
    equi_mdp = make_mdp(mdp = 'equi')
    system(f'{args.gmx} grompp -f {equi_mdp} -c mini.gro -r mini.gro -p topol.top -o equi.tpr -maxwarn 10 {null}')
    deviceID = gpu_manager()
    system(f'{args.mdrun} -deffnm equi -nt {args.numthread} -gpu_id {deviceID}')

#------------- MD
    MD_mdp = make_mdp(mdp = 'MD', ns = args.nanoseconds, dt = args.step)
    system(f'{args.gmx} grompp -f {MD_mdp} -c equi.gro -p topol.top -o MD.tpr -maxwarn 10 {null}')
    deviceID = gpu_manager()
    system(f'{args.mdrun} -v -deffnm MD -nt {args.numthread} -gpu_id {deviceID}')

#------------- ANALISYS
    system(f'{args.gmx} editconf -f MD.tpr -o MD.pdb {null}')
    system(f'echo 0 | {args.gmx} trjconv -f MD.xtc -s MD.tpr -o MD_pbc.xtc -pbc mol')
    system(f"""{args.gmx} mindist -f MD_pbc.xtc -s MD.tpr -d 0.45 -n index.ndx -on contacts.xvg << EOF
    Protein
    13
    EOF""")
    system(f"""{args.gmx} rms -f MD_pbc.xtc -s MD.tpr -n index.ndx << EOF
    24
    13
    EOF""")
    time, all_contacts, all_contacts_fig = plot_xvg('contacts.xvg', 'Number of Contacts',  'contacts.png' ,'Time', 'Concats')
    time, rmsd, _ = plot_xvg('rmsd.xvg', 'RMSD', 'rmsd.png', 'Time', 'RMSD (A)')
    status, contacts_mean = detachmet(all_contacts)
    fig_contacts = contacts(xtc='MD_pbc.xtc', step = args.nanoseconds / 50, ligand = args.ligID)
    result = [filename, status, contacts_mean, all_contacts_fig, np.mean(rmsd) * 10, fig_contacts, smile, platform.node()]
    with open(args.report, 'a') as Report:
        Report.write('\t'.join(map(str, result)) + '\n')
    lines = open(args.report, 'r').readlines()

    if len(lines) % int(args.lignum) == 0:
        df = pd.read_table(args.report, names= columns_name)
        df = df.sort_values(by=['contacts_average', 'status'], ascending = False)
        df = df.head(100)
        df['Molecule'] = df['smiles'].apply(MolFromSmiles)
        df.to_html('../report.html', escape = False)

        send_mail(
                destination='ettore.locascio@unicatt.it',
                subject = f"Report of Ligand-PALS1 {len(lines)}/{len(ligList)}",
                content= f"""This is an auto-generated email from HTMD from {platform.node()}.
                {len(lines)}/{len(ligList)} Processed!
                """,
                attachment = '../report.html',
                )
    else:
        pass
    chdir('../')

if __name__=='__main__':
    p = Pool(int(args.pool))
    p.map_async(main, ligList, chunksize=1).get()
    p.close()
