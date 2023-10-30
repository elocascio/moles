import MDAnalysis
import os
from lists import *
import pickle
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# from scipy.stats import ks_2samp
import argparse


def DNAze_pdb(pdb):
    """If plumed doesn't recognize your DNA, DNAze it!
    replace CYT, ADE, GUA, THY into DC, DA, DG, DT. Ekein"""
    with open(pdb, "rt") as ref:
        data = ref.read()
        data = (
            data.replace("CYT", " DC")
            .replace("ADE", " DA")
            .replace("GUA", " DG")
            .replace("THY", " DT")
        )
    with open(pdb, "wt") as ref:
        ref.write(data)
    print(f"{pdb} resolved")


parser = argparse.ArgumentParser()
parser.add_argument(
    "-tpr", 
    type=str, 
    help="tpr file to use as reference", 
    default="./MD.tpr"
)
parser.add_argument(
    "-xtc", 
    type=str, 
    help="xtc file to use as reference", 
    default="./MD.xtc"
)
parser.add_argument(
    "-ref", 
    type=str, 
    help="path of matrix referement"
)
parser.add_argument(
    "-pdb", 
    type=str, 
    help="selection group A", 
    default="./MD.pdb"
)
args = parser.parse_args()

U = []
df_list = []

# df_ref = '/home/ettore/Documenti/2021/TCF4/TCF4_wt/torsions.df','/home/catt2/Documenti/ETT/torsions.df'
if not os.path.isfile("MD.pdb"):
    os.system(f"gmx editconf -f {args.tpr} -o MD.pdb")

DNAze_pdb(args.pdb)

u = MDAnalysis.Universe(args.pdb)
print(f"PDB file given contains following chains: {u.segments.segids}")
dat = open(f"{os.path.dirname(args.pdb)}/torsion.dat", "w")
colvar = os.path.dirname(args.pdb) + "/colvar_torsion"
dat.write(f"MOLINFO STRUCTURE={args.pdb}\n")
arg_names = []
for chain in u.segments.segids[:-1]:
    n = chain_to_num[chain]
    seg = u.segments[n].residues
    res = seg.resnames[1:-1]
    resnum = seg.resnums[1:-1]
    topol = list(zip(res, resnum))
    for res, resnum in topol:
        arg_name = []
        for angle in X[res]:
            arg_name.append(f"{angle}_{chain}_{resnum}_{res}")
            arg_names.append(f"{angle}_{chain}_{resnum}_{res}")
        zippato = list(zip(X[res], arg_name))
        for angle, arg_name in zippato:
            dat.write(f"""{arg_name}: TORSION ATOMS=@{angle}-{chain}_{resnum}\n""")

dat.write(f"""\nPRINT ARG={','.join(arg_names)} FILE={colvar} STRIDE=10""")
dat.close()
print(dat.name, "written...")
print("PLUMED ANALIZZA", args.xtc, dat.name)
os.system(f"plumed driver --mf_xtc {args.xtc} --plumed {dat.name}")
colvar_file = open(colvar)
data = []
for i, line in enumerate(colvar_file):
    if i == 0:
        fields = line.split()
    elif line.startswith("#") or line.startswith("@"):
        pass
    else:
        x = line.split()
        data.append(x)
colvar_file.close()

df_spec = pd.DataFrame(data, columns=fields[2:])
df_ref = pd.read_pickle(args.ref)

df_list.append(df_ref)
df_list.append(df_spec)

torsion_folder = os.path.join(os.path.dirname(args.pdb), "torsion_plot")
os.makedirs(torsion_folder, exist_ok=True)
os.chdir(torsion_folder)

for chain in u.segments.segids[:-1]:
    n = chain_to_num[chain]
    seg = u.segments[n].residues
    res = seg.resnames[1:-1]
    resnum = seg.resnums[1:-1]
    topol = list(zip(res, resnum))

    for res, resnum in topol:
        try:
            print("ANALIZZANDO: {}_{}_{}".format(chain, resnum, res))
            ref = df_list[0].filter(like="{}_{}_{}".format(chain, resnum, res))
            spec = df_list[1].filter(like="{}_{}_{}".format(chain, resnum, res))
            concat = pd.concat([ref, spec])
            print("TSNE")
            tsne = TSNE(init="pca").fit_transform(concat.values)
            print("plotting {}_{}_{}".format(chain, resnum, res))
            plt.xlabel("PC1")
            plt.ylabel("PC2")
            plt.xticks(range(-50, 51, 10))
            plt.yticks(range(-50, 51, 10))
            plt.title("{}_{}_{}".format(chain, resnum, res))
            plt.scatter(
                tsne[: ref.shape[0], 0],
                tsne[: ref.shape[0], 1],
                c="r",
                marker="+",
                alpha=0.5,
            )
            plt.scatter(
                tsne[ref.shape[0] :, 0],
                tsne[ref.shape[0] :, 1],
                c="blue",
                marker="+",
                alpha=0.5,
            )
            plt.scatter(tsne[0, 0], tsne[0, 1], c="green", s=50, marker="o", alpha=0.9)
            plt.savefig(f"{chain}_{resnum}_{res}.png", format="png")
            plt.close()
        except:
            print(f"MUTATION: {chain}_{resnum}_{res} is mutated!?")
            pass
