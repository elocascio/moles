###########################################
# Runs protein-protein interaction.
#
# Outputs bar plot of the protein-protein interaction
# and data table. Requires topology file for reference (.tpr,.pdb,.psf)
# and a trajectory (.trr,.xtc,.dcd) and selection string (es. "segid PROB and resnum 1-11").
#
# USAGE: python protein_coord.py -h
#
# ettore.locascio@unicatt.it -> version 2023
###########################################
#
# see:
# Giacon, N; Lo Cascio, E et al, Monomeric and dimeric states of human ZO1-PDZ2 are functional partners of the SARS-CoV-2 E protein, Comput Struct Biotechnol J, 2023.
# DOI:https://doi.org/10.1016/j.csbj.2023.05.027
###########################################

def switch(r, r_0=6, a=6, b=12):
    return (1 - (r / r_0) ** a) / (1 - (r / r_0) ** b)


def switch_functions(r, r0, r_0=6, a=6, b=12):
    r = np.asarray(r)
    return np.sum((1 - (r / r_0) ** a) / (1 - (r / r_0) ** b))


import MDAnalysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import argparse


def protein_contacts(topol="MD.pdb", trj="MD.xtc", step=10, ligand="segid PROB"):

    Hydrophobic = pd.DataFrame()
    Hydrophobic_time = pd.DataFrame()
    PiStacking_T = pd.DataFrame()
    PiStacking_P = pd.DataFrame()
    PiCation = pd.DataFrame()
    SaltBridge = pd.DataFrame()
    SaltBridge_time = pd.DataFrame()
    Polar = pd.DataFrame()
    Polar_time = pd.DataFrame()
    HBond = []
    coord_dfs = []

    u = MDAnalysis.Universe(topol, trj)

    for ts in u.trajectory[:step]:

        ########### Hydrophobic ###############
        df_H = pd.DataFrame()
        lig_c = u.select_atoms(f"{ligand} and (name C*)")
        prot_c = u.select_atoms(f"(around 7 {ligand}) and name C*")
        gen = itertools.cycle(range(len(lig_c)))
        prot_list = []
        for res, num in list(zip(prot_c.resnames, prot_c.resnums)):
            prot_list.append(str(num) + str(res))
        for position in lig_c.positions:
            scores = switch(np.linalg.norm(prot_c.positions - position, axis=1))
            df_H[next(gen)] = scores
        df_H = pd.DataFrame((prot_list, df_H.sum(axis=1))).T
        df_H[1] = df_H[1].astype(float)
        df_H = df_H.groupby(0).mean()
        Hydrophobic = pd.concat([Hydrophobic, df_H])
        Hydrophobic_time = pd.concat([Hydrophobic_time, df_H], axis=1)
        ########### Hydrophobic ###############
        ############## Salt bridge ##############
        df_Sn = pd.DataFrame()
        lig_pos = u.select_atoms(f"{ligand} and (name CZ NZ)")
        lig_neg = u.select_atoms(
            f"{ligand} and (resname GLU ASP) and name OT1 OT2 OD1 OD2 OE1 OE2"
        )
        lig_asn = u.select_atoms(f"{ligand} and (resname ASN)")
        lig_neg = lig_neg - lig_asn
        lig_c_term = u.select_atoms(
            f"{ligand} and (resnum {lig_c.resnums[-1]}) and (name OT*)"
        )
        lig_neg = lig_neg + lig_c_term
        prot_pos = u.select_atoms(
            f"(around 7 {ligand}) and (resname ARG LYS) and (name CZ NZ)"
        )
        prot_neg = u.select_atoms(f"(around 7 {ligand}) and (name OD1 OD2 OE1 OE2)")
        # pos vs neg #
        gen = itertools.cycle(range(len(lig_neg)))
        prot_pos_list = []
        for res, num in list(zip(prot_pos.resnames, prot_pos.resnums)):
            prot_pos_list.append(str(num) + str(res))
        for position in lig_neg.positions:
            scores = switch(
                np.linalg.norm(prot_pos.positions - position, axis=1), r_0=5
            )
            df_Sn[next(gen)] = scores
        df_Sn = pd.DataFrame((prot_pos_list, df_Sn.sum(axis=1))).T
        df_Sn[1] = df_Sn[1].astype(float)
        df_Sn = df_Sn.groupby(0).mean()
        # neg vs pos #
        df_Sp = pd.DataFrame()
        gen = itertools.cycle(range(len(lig_pos)))
        prot_neg_list = []
        for res, num in list(zip(prot_neg.resnames, prot_neg.resnums)):
            prot_neg_list.append(str(num) + str(res))
        for position in lig_pos.positions:
            scores = switch(
                np.linalg.norm(prot_neg.positions - position, axis=1), r_0=5
            )
            df_Sp[next(gen)] = pd.Series(scores)
        df_Sp = pd.DataFrame((prot_neg_list, df_Sp.sum(axis=1))).T
        df_Sp[1] = df_Sp[1].astype(float)
        df_Sp = df_Sp.groupby(0).mean()
        df_S = pd.DataFrame()
        df_S = pd.concat([df_Sn, df_Sp])
        SaltBridge = pd.concat([SaltBridge, df_S])
        SaltBridge_time = pd.concat([SaltBridge_time, df_S], axis=1)
        ############## Salt bridge ##############
        ############## Polar ###############
        pos = []
        df_pol = pd.DataFrame()
        lig_back_polar = u.select_atoms(f"{ligand} and (backbone) and (name O* N*)")
        lig_polar = u.select_atoms(
            f"{ligand} and (resname THR SER GLN ASN CYS MET HSD TRP ) and (name O* N* S*)"
        )
        lig_polar = lig_polar + lig_back_polar
        prot_pos = u.select_atoms(
            f"(around 7 {ligand}) and (resname ARG LYS) and (name CZ NZ)"
        )
        gen = itertools.cycle(pos)
        for res, num in list(zip(prot_pos.resnames, prot_pos.resnums)):
            pos.append(str(num) + str(res))
        for position in lig_polar.positions:
            scores = switch(np.linalg.norm(prot_pos.positions - position, axis=1))
            df_pol[next(gen)] = pd.Series(scores)
        df_pol = pd.DataFrame((pos, df_pol.sum(axis=1))).T
        df_pol[1] = df_pol[1].astype(float)
        df_pol = df_pol.groupby(0).mean()

        neg = []
        df_pol2 = pd.DataFrame()
        lig_polar_pos = u.select_atoms(
            f"{ligand} and (resname ARG LYS) and (name CZ NZ)"
        )
        prot_back = u.select_atoms(
            f"(around 7 {ligand}) and (backbone) and (name O* N*)"
        )
        prot_polar = u.select_atoms(
            f"(around 7 {ligand}) and (resname THR SER GLN ASN CYS MET HSD TRP ) and (name O* N* S*)"
        )
        prot_polar = prot_polar + prot_back
        gen = itertools.cycle(neg)
        for res, num in list(zip(prot_polar.resnames, prot_polar.resnums)):
            neg.append(str(num) + str(res))
        for position in lig_polar_pos.positions:
            scores = switch(np.linalg.norm(prot_polar.positions - position, axis=1))
            df_pol2[next(gen)] = pd.Series(scores)
        df_pol2 = pd.DataFrame((neg, df_pol2.sum(axis=1))).T
        df_pol2[1] = df_pol2[1].astype(float)
        df_pol2 = df_pol2.groupby(0).mean()
        df_Pol = pd.DataFrame()
        df_Pol = pd.concat([df_pol, df_pol2])
        Polar = pd.concat([Polar, df_Pol])
        # Polar_time = pd.concat([Polar, df_Pol], axis = 1)
    ################# Polar #####################
    ############# H-bond #####################
    hbondsa = HBA(universe=u)
    hbondsa.hydrogens_sel = f"(around 7 {ligand}) and (name H* and bonded name O* N*)"
    hbondsa.acceptors_sel = f"{ligand} and (name O* N*)"
    hbondsa.run(step=step)
    resulta = hbondsa.count_by_ids()

    hbondsb = HBA(universe=u)
    hbondsb.hydrogens_sel = f"{ligand} and (name H* and bonded name O* N*)"
    hbondsb.acceptors_sel = f"(around 7 {ligand}) and (name O* N*)"
    hbondsb.run(step=step)
    resultb = hbondsb.count_by_ids()

    for val in resulta:
        HBond.append(
            [str(u.atoms[val[1]].resnum) + str(u.atoms[val[1]].resname), val[3]]
        )
    for val in resultb:
        HBond.append(
            [str(u.atoms[val[2]].resnum) + str(u.atoms[val[2]].resname), val[3]]
        )
    h = pd.DataFrame(HBond)
    h.reset_index()

    ############ FINAL ################
    for df in [Hydrophobic, SaltBridge, df_Pol]:
        df = df.reset_index()
        df.columns = ["residue", "val"]
        df = df.groupby(df["residue"]).aggregate({"val": "sum"})
        df["val"] = df["val"] / df["val"].max()
        coord_dfs.append(df)
    h.columns = ["residue", "val"]
    df = h.groupby(h["residue"]).aggregate({"val": "sum"})
    df["val"] = df["val"] / df["val"].max()
    coord_dfs.append(df)

    df_all = pd.concat(coord_dfs, axis=1)
    df_all = df_all.fillna(0)
    df_all = df_all[(df_all.T > 0.2).any()]
    df_all = df_all.sort_values(by=["residue"])
    df_all.columns = ["Hydrophobic", "SaltBridge", "Polar", "H-bond"]
    df_all.to_csv(f"coord_{ligand}.csv")
    ax = df_all.plot.bar(
        stacked=True, color=["purple", "mediumblue", "lightgreen", "aqua"]
    )
    ax.set_ylim(top=4)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="best")
    ax.set_ylabel("coordination")
    (ax.figure).savefig(f"coord_{ligand}.png", format="png", bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-select", type=str, help="select ligand usign MDAnalsysis selection"
    )
    parser.add_argument(
        "-topol", type=str, help="file pdb,psf,tpr --- default MD.pdb", default="MD.pdb"
    )
    parser.add_argument(
        "-trj", type=str, help="file xtc,trr,dcd --- default MD.xtc", default="MD.xtc"
    )
    parser.add_argument(
        "-step", type=int, help="value of step --- default 10", default="10"
    )
    args = parser.parse_args()
    protein_contacts(topol=args.topol, trj=args.trj, step=args.step, ligand=args.lig)
