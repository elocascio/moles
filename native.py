import MDAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts
import numpy as np
import argparse
from Misc.moles import init
import seaborn as sns
from os import system, listdir, chdir
import pickle as pkl
from os.path import isfile, isdir

init()
parser = argparse.ArgumentParser()
parser.add_argument("-A", type=str)
parser.add_argument("-B", type=str)
parser.add_argument("-method", type=str, choices=['switch_function', 'hard'], )
parser.add_argument("-ref", type=str, help="path of matrix referement" )
parser.add_argument("-sub", action='store_true', help="Calcualte Delta Matrix")
parser.add_argument("-inferno", action='store_true', help="Multidir")
parser.add_argument("-sel_A", type=str, help="selection group A", default=0)
parser.add_argument("-sel_B", type=str, help="selection group B", default=0)
parser.add_argument("-step", type=int, help="step every n frame", default=10)

args = parser.parse_args()



def switch_functions(r, r0, r_0 = 6, a = 6, b = 12):
    r = np.asarray(r)
    return np.sum((1-(r/r_0)**a)/(1-(r/r_0)**b))

if args.inferno == False:
    if not isfile('MD.pdb'):
        system('gmx editconf -f MD.tpr -o MD.pdb')

    u = MDAnalysis.Universe('MD.pdb','MD.xtc')

    if args.sel_A == 0 and args.sel_B == 0:
        target_A_str = f'(segid {args.A}) and not (backbone)' ; target_A = u.select_atoms(target_A_str)
        target_B_str = f'(segid {args.B}) and not (backbone)'; target_B = u.select_atoms(target_B_str)
    else:
        target_A_str = args.sel_A ; target_A = u.select_atoms(target_A_str)
        target_B_str = args.sel_B; target_B = u.select_atoms(target_B_str)
    # CONTACTS
    times = []

    for at in target_B.residues:
        print(at, at.segid)
        prot_str = f"(resid {at.resid}) and (segid {at.segid})"
        coord_P = contacts.Contacts(u,
                                    select = (target_A_str, prot_str) ,
                                    refgroup = (target_A, u.select_atoms(prot_str)) ,
                                    radius = 5 ,
                                    method= 'soft_cut' ,
                                   # kwargs={'r_0':5.5, 'a':6, 'b':12}
                                    ).run(step=args.step)
    #    P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
        times.append(coord_P.timeseries[:,1])

    with open(f'{args.sel_B}_with_{args.sel_A}.pkl', 'wb') as f: pkl.dump(times, f)

    plt.figure(figsize=(30,20))
    ax = sns.heatmap(times, yticklabels=list(zip(target_B.residues.resnums, target_B.residues.resnames)), cmap="Blues")
    ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize = 36)
    plt.savefig(f'contacts_of_{args.sel_B}_with_{args.sel_A}.png', format = 'png')
    plt.close()

    if args.sub:
        matrix = pd.DataFrame(times); matrix = matrix.fillna(0)  
        rmatrix = pkl.load(open(args.ref, 'rb')); rmatrix = pd.DataFrame(rmatrix); rmatrix = rmatrix.fillna(0)
        dmatrix = matrix.subtract(rmatrix)

        plt.figure(figsize=(30,20))
        plt.title(f"delta_{args.sel_A} vs {args.sel_B}")
        ax = sns.heatmap(dmatrix, yticklabels=list(zip(target_B.residues.resnums, target_B.residues.resnames)), cmap="RdBu")
        plt.savefig(f"delta_matrix_{args.B}_with_{args.A}")

#-------------------------------------------------------------------------------
if args.inferno:

    for x in listdir('.'):
        if isdir(x) and isfile(f'{x}/MD.tpr'):

            if not isfile(f'{x}/MD.pdb'):
                system(f'gmx editconf -f {x}/MD.tpr -o {x}/MD.pdb')
            if isfile(f'{x}/contacts_of_{args.sel_B}_with_{args.sel_A}2.png'):
                continue
            try:
                u = MDAnalysis.Universe(f'{x}/MD.pdb',f'{x}/MD.xtc')

                target_A_str = f'{args.sel_A}' ; target_A = u.select_atoms(target_A_str)
                target_B_str = f'{args.sel_B}'; target_B = u.select_atoms(target_B_str)

                # CONTACTS
                times = []
                print(f'Analysing {x}')
                for at in target_B.residues:
                    prot_str = f"(resid {at.resid}) and (segid {at.segid})"
                    coord_P = contacts.Contacts(u,
                                                select = (target_A_str, prot_str) ,
                                                refgroup = (target_A, u.select_atoms(prot_str)) ,
                                                radius = 5 ,
                                                method= 'soft_cut' ,
                                               #kwargs={'r_0':5.5, 'a':6, 'b':12}
                                                ).run(step=10)
                    #P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
                    times.append(coord_P.timeseries[:,1])

                with open(f'{x}/{args.sel_B}_with_{args.sel_A}.pkl', 'wb') as f: pkl.dump(times, f)

                plt.figure(figsize=(30,20))
                ax = sns.heatmap(times, yticklabels=list(zip(target_B.residues.resnums, target_B.residues.resnames)), cmap="RdBu")
                ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize = 36)                
                plt.savefig(f'{x}/contacts_of_{args.sel_B}_with_{args.sel_A}.png', format = 'png')
                plt.close()

                if args.sub:
                    matrix = pd.DataFrame(times); matrix = matrix.fillna(0)  
                    rmatrix = pkl.load(open(args.ref, 'rb')); rmatrix = pd.DataFrame(rmatrix); rmatrix = rmatrix.fillna(0)
                    dmatrix = matrix.subtract(rmatrix)

                plt.figure(figsize=(30,20))
                ax = sns.heatmap(dmatrix, yticklabels=list(zip(target_B.residues.resnums, target_B.residues.resnames)), cmap="RdBu")
                ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize = 36)
                plt.savefig(f"{x}/delta_matrix_{args.sel_B}_with_{args.sel_A}")
            except:
                print('something wrong')