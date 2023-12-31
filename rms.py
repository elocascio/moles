from tkinter import N
import numpy as np
from os import system as sys
import matplotlib.pyplot as plt
import argparse
from utils import clean_xvg
import MDAnalysis as mda
from Misc.moles import init
from biopandas.pdb import PandasPdb as ppdb
import matplotlib as mpl

parser = argparse.ArgumentParser()
parser.add_argument(
    "-tpr",
    type=str,
    help="file tpr gromacs --- default MD.tpr",
    default="MD.tpr",
    required=True,
)
parser.add_argument(
    "-trr",
    type=str,
    help="file trr gromacs --- default MD.trr",
    default="MD.trr",
    required=True,
)
parser.add_argument(
    "-clust", 
    action="store_true", 
    help="Calculate Clusters"
)
parser.add_argument(
    "-method",
    type=str,
    help="Method for cluster determination: linkage, jarvis-patrick, monte-carlo, diagonalization, gromos --- default gromos",
    default="gromos",
)
parser.add_argument(
    "-cutoff", 
    type=float, 
    help="n of cutoff --- default 0.1", 
    default=0.1
)
parser.add_argument(
    "-fit", 
    type=str, 
    help="select group charmm-like syntax for fit-least square group"
)
parser.add_argument(
    "-ogroup", 
    type=str, 
    help="select group charmm-like syntax for group output"
)
parser.add_argument(
    "-skip", 
    type=str, 
    help="skip frame", 
    default=10
)
parser.add_argument(
    "-rmsd", 
    action="store_true", 
    help="Calculate Clusters"
)
parser.add_argument(
    "-stride", 
    type=str, 
    help="frequency", 
    default=10
)
parser.add_argument(
    "-sele",
    type=str,
    help='select group charmm-like syntax for superimposition ex. "name CA or resname LIG"',
)
parser.add_argument(
    "-resname",
    type=str,
    help="resname for beta factor ex. LIG. Leave empty to calculate RMSD of all thing.",
)
args = parser.parse_args()

init()

if args.clust:
    u = mda.Universe(args.tpr, args.trr)
    fit = u.select_atoms(args.fit)
    fit_str = "_".join((args.fit).split())  # protein, name CA, backbone, name UNK
    ogroup = u.select_atoms(args.ogroup)
    ogroup_str = "_".join(
        (args.ogroup).split()
    )  # takes selection string and join with "_"

    with mda.selections.gromacs.SelectionWriter(f"clust.ndx", mode="w") as ndx:
        ndx.write(fit, name=fit_str)
        ndx.write(ogroup, name=ogroup_str)

    sys(
        f"""gmx cluster -f {args.trr} -s {args.tpr} -n clust.ndx -sz cluster_size_{fit_str}.xvg -clid cluster_id_{fit_str}.xvg -cl cluster_repr_{fit_str}.pdb -method {args.method} -cutoff {args.cutoff} -skip {args.skip} << EOF
    0
    1
    EOF"""
    )

    cmap = mpl.colors.ListedColormap(
        ["blue", "orange", "green", "red", "pink", "purple", "grey", "cyan"]
    )

    bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    n_clust, freq = clean_xvg(f"cluster_size_{fit_str}.xvg")
    colors = ["blue", "orange", "green", "red", "pink", "purple", "grey", "cyan"]
    if len(freq) > 7:
        n = 7
        plt.title(f"CLUSTER {fit_str}")
        others = np.sum(freq[7:])
        freq_new = freq[:7]
        print(freq_new)
        freq_new = np.append(freq_new, others)
        print(freq_new)
        plt.pie(freq_new, colors=colors[:8])
        labels = []
        for ind, q in list(zip(range(len(freq_new)), freq_new)):
            if ind < 7:
                labels.append(f"{ind + 1}: {np.round(q/freq_new.sum()*100, 2)}")
            else:
                labels.append(f"OTHERS: {np.round(q/freq_new.sum()*100, 2)}")
        plt.legend(labels[:10], title="CLUSTER %", bbox_to_anchor=(1, 0.7))
        plt.tight_layout()
        plt.savefig(
            f"clust_{fit_str}_{args.cutoff}_{args.method}.png", format="png", dpi=900
        )
        plt.close()
    else:
        n = len(freq)
        labels = []
        for ind, q in list(zip(range(n), freq)):
            labels.append(f"{ind + 1}: {np.round(q/freq.sum()*100, 2)}")
        plt.title(f"CLUSTER {fit_str}")
        plt.pie(freq[:8], colors=colors[:n])
        plt.legend(labels[:10], title="CLUSTER %", bbox_to_anchor=(1, 0.7))
        plt.tight_layout()
        plt.savefig(
            f"clust_{fit_str}_{args.cutoff}_{args.method}.png", format="png", dpi=900
        )
        plt.close()

    time, clust = clean_xvg(f"cluster_id_{fit_str}.xvg")
    fig, ax = plt.subplots()
    ax.imshow([clust, clust, clust], cmap=cmap, norm=norm)
    ax.set_title("Cluster")
    ax.set_xticks([0, len(clust) / 2, len(clust)])
    ax.set_xticklabels([0, round(time[-1] / 2000), round(time[-1] / 1000)])
    plt.savefig(f"cluster_map_{fit_str}.png", format="png", dpi=900)
    plt.close()

if args.rmsd:
    u = mda.Universe(args.tpr, args.trr)
    selection = u.select_atoms(args.sele)  # protein, name CA, backbone, name UNK
    sele = "_".join((args.sele).split())  # takes selection string and join with "_"
    with mda.selections.gromacs.SelectionWriter(f"{sele}.ndx", mode="w") as ndx:
        ndx.write(selection, name=sele)
    sys(f"gmx editconf -f {args.tpr} -o ref_{sele}.pdb -n {sele}.ndx")

    if args.resname:
        file = ppdb().read_pdb(f"ref_{sele}.pdb")
        file.df["ATOM"]["b_factor"] = np.where(
            (file.df["ATOM"]["residue_name"] == (args.resname)[:3]),
            1,
            file.df["ATOM"]["b_factor"],
        )
        file.to_pdb(path=f"ref_{sele}.pdb")

    with open("rmsd.dat", "w") as dat:
        dat.write(
            f"""
UNITS LENGTH=A
RMSD REFERENCE=ref_{sele}.pdb TYPE=OPTIMAL
PRINT ARG=* FILE=rmsd_{'_'.join((args.resname).split())} STRIDE={args.stride}"""
        )

    sys(f"plumed driver --mf_{args.trr[-3:]} {args.trr} --plumed rmsd.dat")

    time, rmsd = clean_xvg(f"rmsd_{'_'.join((args.resname).split())}")

    plt.title(f"RMSD {args.resname}")
    plt.plot(time, rmsd)
    plt.xlabel("Time(ns)")
    plt.ylabel("RMSD (A°)")
    # plt.legend(labels,title="CLUSTER %",bbox_to_anchor=(1, .7))
    plt.tight_layout()
    plt.savefig(f"RMSD_{args.resname}.png", format="png", dpi=900)
