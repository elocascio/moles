
import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f",type=str, help="file tpr gromacs --- default MD.tpr", default= 'MD.tpr', required=True)
parser.add_argument("-sele", type=str, help="select group charmm-like syntax for superimposition ex. \"name CA or resname LIG\"")
args = parser.parse_args()


u = mda.Universe(args.f)
selection = u.select_atoms(args.sele) # protein, name CA, backbone, name UNK
sele = '_'.join((args.sele).split()) # takes selection string and join with "_"
with mda.selections.gromacs.SelectionWriter(f'{sele}.ndx', mode = 'w') as ndx:
    ndx.write(selection, name = sele)