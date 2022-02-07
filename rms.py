import numpy as np 
from os import system as sys
import matplotlib.pyplot as plt
import argparse
from utils import clean_xvg
import MDAnalysis as mda

parser = argparse.ArgumentParser()
parser.add_argument("-clust", action="store_true", help="Calculate Clusters")
parser.add_argument("-rmsd", action="store_true", help="Calculate Clusters")
parser.add_argument("-trr",type=str, help="file trr gromacs --- default MD.trr", default='MD.trr')
parser.add_argument("-tpr",type=str, help="file tpr gromacs --- default MD.tpr", default= 'MD.tpr')
parser.add_argument("-ndx",type=str, help="file index gromacs --- default index.ndx", default='index.ndx')
parser.add_argument("-method",type=str, help="Method for cluster determination: linkage, jarvis-patrick, monte-carlo, diagonalization, gromos --- default gromos",default='gromos')
parser.add_argument("-cutoff",type=float, help="n of cutoff --- default 0.1", default=0.1,)
parser.add_argument("-fit", type=str, help="fit-least square group")
parser.add_argument("-ogroup", type=str, help="ogroup-group output")
parser.add_argument("-stride", type=str, help="frequency")
parser.add_argument("-sele", type=str, help="select group charmm-like syntax") 
 
args = parser.parse_args()


if args.clust:
    sys(f'''gmx cluster -f {args.trr} -s {args.tpr} -n {args.ndx} -sz -cl -method {args.method} -cutoff {args.cutoff} << EOF
    {args.fit}
    {args.ogroup}
    EOF''')

    n_clust,freq = clean_xvg('clust-size.xvg')

    fig = plt.figure(figsize=(5,7))
    plt.title('CLUSTER')
    plt.pie(freq)
    labels= np.round(freq/freq.sum()*100, 2)
    plt.legend(labels,title="CLUSTER %",bbox_to_anchor=(1, .7))
    plt.tight_layout()
    fig.savefig(f'clust_{args.cutoff}_{args.method}.png',format='png', dpi=1600)

if args.rmsd:
    u = mda.Universe(args.tpr, args.trr)
    selection = u.select_atoms(args.sele) # protein, name CA, backbone, name UNK
    selection.write(f'reference_rmsd.pdb')
    with open('rmsd.dat', 'w') as dat:
        dat.write(f"""
RMSD REFERENCE=reference_rmsd.pdb TYPE=OPTIMAL
PRINT ARG=* FILE=rmsd STRIDE={args.stride}""")

    sys(f'plumed driver --mf_{args.trr[-3:]} {args.trr} --plumed rmsd.dat')

    time, rmsd = clean_xvg('rmsd')

    fig = plt.figure(figsize=(11,7))
    plt.title(f'RMSD {args.sele}')
    plt.plot(time, rmsd)
    plt.xlabel('Time(ns)')
    plt.ylabel('RMSD (A°)')
    #plt.legend(labels,title="CLUSTER %",bbox_to_anchor=(1, .7))
    plt.tight_layout()
    fig.savefig(f'RMSD_{args.sele}.png',format='png', dpi=1600)
    