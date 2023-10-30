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
import sys

init()
parser = argparse.ArgumentParser()
parser.add_argument(
    "-tpr", type=str, help="tpr file to use as reference", default="MD.tpr"
)
parser.add_argument(
    "-xtc", type=str, help="xtc file to use as reference", default="MD.xtc"
)

parser.add_argument("-method", type=str, choices=["switch_function", "hard"])
parser.add_argument("-ref", type=str, help="path of matrix referement")
parser.add_argument("-sub", action="store_true", help="Calcualte Delta Matrix")
parser.add_argument("-multi", action="store_true", help="Multidir")
parser.add_argument("-show", action="store_true", help="show")
parser.add_argument("-sel_A", type=str, help="selection group A", default=0)
parser.add_argument("-sel_B", type=str, help="selection group B", default=0)
parser.add_argument("-step", type=int, help="step every n frame", default=10)
parser.add_argument(
    "-absolute", action="store_true", help="absolute score using switching function"
)
parser.add_argument("-div", action="store_true", help="divide 2 absolute-scored matrix")

args = parser.parse_args()


def switch_functions(r, r0, r_0=5, a=6, b=12):
    r = np.asarray(r)
    return np.sum((1 - (r / r_0) ** a) / (1 - (r / r_0) ** b))


if args.multi == False:
    if not isfile("MD.pdb"):
        system(f"gmx editconf -f {args.tpr} -o MD.pdb")

    u = MDAnalysis.Universe("MD.pdb", args.xtc)

    target_A_str = f"{args.sel_A}"
    target_A = u.select_atoms(target_A_str)
    target_B_str = f"{args.sel_B}"
    target_B = u.select_atoms(target_B_str)

    # CONTACTS
    times = []
    if args.absolute == False:
        for at in target_B.residues:
            prot_str = f"(resid {at.resid} and (segid {at.segid}))"
            coord_P = contacts.Contacts(
                u,
                select=(target_A_str, prot_str),
                refgroup=(target_A, u.select_atoms(prot_str)),
                radius=5,
                method="soft_cut",
                # kwargs={'r_0':5.0, 'a':6, 'b':12}
            ).run(step=args.step)
            # P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
            times.append(coord_P.timeseries[:, 1])
        with open(f"./{args.sel_B}_with_{args.sel_A}.pkl", "wb") as f:
            pkl.dump(times, f)
        np.savetxt("times", times)
        plt.figure(figsize=(30, 20))
        ax = sns.heatmap(
            times,
            yticklabels=list(
                zip(target_B.residues.resnums, target_B.residues.resnames)
            ),
            cmap="Blues",
            vmin=0,
            vmax=1,
        )
        ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize=36)
        plt.savefig(f"./contacts_of_{args.sel_B}_with_{args.sel_A}.png", format="png")
        plt.close()
        print("matrix plot saved")
    else:
        for at in target_B.residues:
            prot_str = f"(resid {at.resid} and (segid {at.segid}))"
            coord_P = contacts.Contacts(
                u,
                select=(target_A_str, prot_str),
                refgroup=(target_A, u.select_atoms(prot_str)),
                radius=10,
                method=switch_functions,
                kwargs={"r_0": 6, "a": 6, "b": 12},
            ).run(step=args.step)
            # P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
            times.append(coord_P.timeseries[:, 1])
        with open(f"./abs{args.sel_B}_with_{args.sel_A}.pkl", "wb") as f:
            pkl.dump(times, f)
        plt.figure(figsize=(30, 20))
        ax = sns.heatmap(
            times,
            yticklabels=list(
                zip(target_B.residues.resnums, target_B.residues.resnames)
            ),
            cmap="Blues",
            vmin=0,
        )
        ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize=36)
        plt.savefig(
            f"./abs_contacts_of_{args.sel_B}_with_{args.sel_A}.png", format="png"
        )
        plt.close()
        # plt.plot(coord_P.timeseries[:,0], coord_P.timeseries[:,1]); plt.savefig('plot.png', fomrmat = 'png'); plt.close()
        print("matrix plot saved")
    if args.sub:
        matrix = pd.DataFrame(times)
        matrix = matrix.fillna(0)
        rmatrix = pkl.load(open(args.ref, "rb"))
        rmatrix = pd.DataFrame(rmatrix)
        rmatrix = rmatrix.fillna(0)
        dmatrix = matrix.subtract(rmatrix)
        dmatrix["mean"] = dmatrix.mean(axis=1)
        dmatrix["residue"] = list(
            zip(target_B.residues.resnums, target_B.residues.resnames)
        )
        selected = dmatrix[dmatrix["mean"] < -0.3]
        print(selected[["residue", "mean"]])
        plt.figure(figsize=(30, 20))
        plt.title(f"delta_{args.sel_A} vs {args.sel_B}")
        ax = sns.heatmap(
            dmatrix,
            yticklabels=list(
                zip(target_B.residues.resnums, target_B.residues.resnames)
            ),
            cmap="RdBu",
            vmin=-1,
            vmax=1,
        )
        plt.savefig(f"delta_matrix_{args.sel_B}_with_{args.sel_A}")

# -------------------------------------------------------------------------------
if args.multi:

    for x in listdir("."):
        if isdir(x) and isfile(f"{x}/MD.tpr"):

            if not isfile(f"{x}/MD.pdb"):
                system(f"gmx editconf -f {x}/MD.tpr -o {x}/MD.pdb")
            try:
                u = MDAnalysis.Universe(f"{x}/MD.pdb", f"{x}/MD.xtc")

                target_A_str = f"{args.sel_A}"
                target_A = u.select_atoms(target_A_str)
                target_B_str = f"{args.sel_B}"
                target_B = u.select_atoms(target_B_str)

                # CONTACTS
                times = []
                print(f"Analysing {x}")
                if args.absolute == False:
                    for at in target_B.residues:
                        prot_str = f"(resid {at.resid} and (segid {at.segid}))"
                        coord_P = contacts.Contacts(
                            u,
                            select=(target_A_str, prot_str),
                            refgroup=(target_A, u.select_atoms(prot_str)),
                            radius=5,
                            method="soft_cut",
                            # kwargs={'r_0':5.0, 'a':6, 'b':12}
                        ).run(step=args.step)
                        # P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
                        times.append(coord_P.timeseries[:, 1])
                        with open(f"{x}/{args.sel_B}_with_{args.sel_A}.pkl", "wb") as f:
                            pkl.dump(times, f)
                        plt.figure(figsize=(30, 20))
                    ax = sns.heatmap(
                        times,
                        yticklabels=list(
                            zip(target_B.residues.resnums, target_B.residues.resnames)
                        ),
                        cmap="Blues",
                        vmin=0,
                        vmax=1,
                    )
                    ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize=36)
                    plt.savefig(
                        f"{x}/contacts_of_{args.sel_B}_with_{args.sel_A}.png",
                        format="png",
                    )
                    plt.close()
                    print("matrix plot saved")
                else:
                    for at in target_B.residues:
                        prot_str = f"(resid {at.resid} and (segid {at.segid}))"
                        coord_P = contacts.Contacts(
                            u,
                            select=(target_A_str, prot_str),
                            refgroup=(target_A, u.select_atoms(prot_str)),
                            radius=10,
                            method=switch_functions,
                            kwargs={"r_0": 4.5, "a": 6, "b": 12},
                        ).run(step=args.step)
                        # P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
                        times.append(coord_P.timeseries[:, 1])
                    with open(f"{x}/abs{args.sel_B}_with_{args.sel_A}.pkl", "wb") as f:
                        pkl.dump(times, f)

                    plt.figure(figsize=(30, 20))
                    ax = sns.heatmap(
                        times,
                        yticklabels=list(
                            zip(target_B.residues.resnums, target_B.residues.resnames)
                        ),
                        cmap="Blues",
                        vmin=0,
                    )
                    ax.set_title(f"{args.sel_A} vs {args.sel_B}", fontsize=36)
                    plt.savefig(
                        f"{x}/abs_contacts_of_{args.sel_B}_with_{args.sel_A}.png",
                        format="png",
                    )
                    plt.close()
                    print("matrix plot saved")

                if args.div and args.absolute:
                    print("NOW DIV")
                    matrix = pd.DataFrame(times)
                    matrix = matrix.fillna(0)
                    rmatrix = pkl.load(open(args.ref, "rb"))
                    rmatrix = pd.DataFrame(rmatrix)
                    rmatrix = rmatrix.fillna(0)
                    dmatrix = matrix.div(rmatrix)
                    dmatrix = dmatrix.fillna(0)
                    plt.figure(figsize=(30, 20))
                    ax = sns.heatmap(
                        dmatrix,
                        yticklabels=list(
                            zip(target_B.residues.resnums, target_B.residues.resnames)
                        ),
                        cmap="RdBu",
                        vmin=0,
                        vmax=1,
                    )
                    ax.set_title(f"{args.sel_B} vs {args.sel_A}", fontsize=36)
                    plt.savefig(
                        f"{x}/div_abs_matrix_{args.sel_B}_with_{args.sel_A}.png",
                        format="png",
                    )
                    plt.close()
                    dmatrix["mean"] = dmatrix.mean(axis=1)
                    dmatrix["residue"] = list(
                        zip(target_B.residues.resnums, target_B.residues.resnames)
                    )
                    selected = dmatrix[dmatrix["mean"] < 0.6]
                    print(selected[["residue", "mean"]])
                elif args.sub:
                    print("NOW SUB")
                    matrix = pd.DataFrame(times)
                    matrix = matrix.fillna(0)
                    rmatrix = pkl.load(open(args.ref, "rb"))
                    rmatrix = pd.DataFrame(rmatrix)
                    rmatrix = rmatrix.fillna(0)
                    dmatrix = matrix.subtract(rmatrix)
                    plt.figure(figsize=(30, 20))
                    if args.absolute:
                        ax = sns.heatmap(
                            dmatrix,
                            yticklabels=list(
                                zip(
                                    target_B.residues.resnums,
                                    target_B.residues.resnames,
                                )
                            ),
                            cmap="RdBu",
                        )
                    else:
                        ax = sns.heatmap(
                            dmatrix,
                            yticklabels=list(
                                zip(
                                    target_B.residues.resnums,
                                    target_B.residues.resnames,
                                )
                            ),
                            cmap="RdBu",
                            vmin=-1,
                            vmax=1,
                        )
                    ax.set_title(f"{args.sel_B} vs {args.sel_A}", fontsize=36)

                    if args.absolute:
                        plt.savefig(
                            f"{x}/abs_delta_matrix_{args.sel_B}_with_{args.sel_A}.png",
                            format="png",
                        )
                    else:
                        plt.savefig(
                            f"{x}/delta_matrix_{args.sel_B}_with_{args.sel_A}.png",
                            format="png",
                        )
                    plt.close()
                    dmatrix["mean"] = dmatrix.mean(axis=1)
                    dmatrix["residue"] = list(
                        zip(target_B.residues.resnums, target_B.residues.resnames)
                    )
                    selected = dmatrix[dmatrix["mean"] < 0.1]
                    print(selected[["residue", "mean"]])
                else:
                    print("what?")
            except OSError as err:
                print("OS error: {0}".format(err))
            except ValueError:
                print("Could not convert data to an integer.")
            except:
                print("Unexpected error:", sys.exc_info()[0])
