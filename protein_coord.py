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

import MDAnalysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import argparse

def _parse_CLAs() -> argparse.ArgumentParser:
    """
    Parser Command Line Arguments

    Returns:
        argparse.ArgumentParser : parsed command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-select", 
        type=str, 
        help="select ligand usign MDAnalsysis selection"
    )
    parser.add_argument(
        "-topol", 
        type=str, 
        help="file pdb,psf,tpr --- default MD.pdb",
        default="MD.pdb"
    )
    parser.add_argument(
        "-trj", 
        type=str, 
        help="file xtc,trr,dcd --- default MD.xtc", 
        default="MD.xtc"
    )
    parser.add_argument(
        "-step", 
        type=int, 
        help="value of step --- default 10", 
        default="10"
    )
    parser.add_argument(
        "-limit", 
        type=float, 
        help="limit value for interaction histogram --- default .2", 
        default=".2"
    )
    return parser

def switch(r: np.ndarray, r_0: float =6, a: int =6, b: int =12) -> np.ndarray:
    """
    Calculate interaction scores using a rational switching function from PLUMED.

    Args:
        r (np.ndarray): Atom-atom distances.
        r_0 (float): Parameter of the switching function (default is 6).
        a (int): Exponent for the numerator (default is 6).
        b (int): Exponent for the denominator (default is 12).

    Returns:
        np.ndarray: Interaction scores converted from distances using the switching function.
    """
    return (1 - (r / r_0) ** a) / (1 - (r / r_0) ** b)

def calculate_center(points: np.ndarray) -> np.ndarray:
    """
    calculate center from array of points
    
    Args:
        points (np.ndarray): coordinates of each point
    
    Returns:
        np.ndarray: coordinates of the central point
        """
    sum = np.sum(points, axis = 0)
    center = sum / len(points)
    return center

def calculate_plane_normal(p1, p2, p3):
    """
    Calculate the normal vector to a plane given three points.

    Args:
        p1 (np.ndarray): Coordinates of the first point.
        p2 (np.ndarray): Coordinates of the second point.
        p3 (np.ndarray): Coordinates of the third point.

    Returns:
        np.ndarray: Normal vector to the plane.
    """
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)

def calculate_angle_between_planes(plane1_points, plane2_points):
    """
    Calculate the angle in degrees between two planes given three points each.

    Args:
        plane1_points (tuple): Tuple of three points identifying the first plane.
        plane2_points (tuple): Tuple of three points identifying the second plane.

    Returns:
        float: Angle in degrees between the two planes.
    """
    normal1 = calculate_plane_normal(*plane1_points)
    normal2 = calculate_plane_normal(*plane2_points)
    dot_product = np.dot(normal1, normal2)
    angle_radians = np.arccos(dot_product)
    angle_degrees = np.degrees(angle_radians)
    return angle_degrees

def calculate_distance_and_angle_to_plane(point, plane_points):
    """
    Calculate the distance between a point and a plane, and the angle between the plane and the point.

    Args:
        point (np.ndarray): Coordinates of the point.
        plane_points (tuple): Tuple of three points identifying the plane.

    Returns:
        tuple: A tuple containing the distance between the point and the plane, and the angle in degrees.
    """
    # Calculate the normal vector to the plane
    plane_normal = calculate_plane_normal(*plane_points)
    center = calculate_center(*plane_points)
    # Calculate the vector from any point on the plane to the given point
    vector_to_point = point - plane_points[0]

    # Calculate the distance between the point and the plane
    distance_to_plane = np.linalg.norm(center - point)

    # Calculate the angle between the plane and the vector to the point
    angle_to_plane = np.arccos(np.dot(vector_to_point, plane_normal) / (np.linalg.norm(vector_to_point) * np.linalg.norm(plane_normal)))
    angle_to_plane_degrees = np.degrees(angle_to_plane)

    return distance_to_plane, angle_to_plane_degrees

def diction_me(res_list, values, dictionary):
    """
    Create or update a dictionary of interactions using the provided residue list and corresponding values.

    Args:
        res_list (iterable): List of residue identifiers.
        values (iterable): List of corresponding interaction values.
        dictionary (dict): Existing dictionary of interactions to update, or an empty dictionary.

    Returns:
        dict: Updated dictionary of interactions where residue identifiers (keys) are associated with their cumulative values.
    """
    for res, score in zip(res_list, values):
        if res in dictionary:
            dictionary[res] += score
        else:
            dictionary[res] = score

    return dictionary

def calculate_hydrophobic_interaction(u, ligand, hydrophobic_data = {}):
    """
    Calculate hydrophobic interaction scores between a ligand and surrounding protein atoms.

    Args:
        u (MDAnalysis.Universe): MDAnalysis Universe object representing the molecular system.
        ligand (str): Selection string for the ligand of interest.
        hydrophobic_data (dict, optional): Dictionary to store hydrophobic interaction scores (default is an empty dictionary).

    Returns:
        dict: Updated dictionary of hydrophobic interaction scores, where each protein residue is associated with its cumulative score.
    """

    lig_c = u.select_atoms(f"{ligand} and (name C*)")
    prot_c = u.select_atoms(f"(around 7 {ligand}) and name C*")
    prot_c_res = [str(num) + str(res) for res, num in zip(prot_c.resnames, prot_c.resnums)]

    for position in lig_c.positions:
        scores = switch(np.linalg.norm(prot_c.positions - position, axis=1))
    
    hydrophobic_data = diction_me(prot_c_res, scores, hydrophobic_data)

    return hydrophobic_data

def calculate_salt_bridge_interaction(u, ligand, salt_data = {}):
    """
    Calculate salt bridge interaction scores between a ligand and surrounding protein atoms.

    Args:
        u (MDAnalysis.Universe): MDAnalysis Universe object representing the molecular system.
        ligand (str): Selection string for the ligand of interest.
        salt_data (dict, optional): Dictionary to store salt bridge interaction scores (default is an empty dictionary).

    Returns:
        dict: Updated dictionary of salt bridge interaction scores, where each protein residue is associated with its cumulative score.
    """

    lig_pos = u.select_atoms(f"{ligand} and (name CZ NZ)")
    lig_neg = u.select_atoms(f"{ligand} and (name OT1 OT2 OD1 OD2 OE1 OE2) and not (resname ASN GLN)")
    prot_pos = u.select_atoms(f"(around 7 {ligand}) and (resname ARG LYS) and (name CZ NZ)")
    prot_neg = u.select_atoms(f"(around 7 {ligand}) and (name OD1 OD2 OE1 OE2) and not (resname ASN GLN)")
    
    # protein pos vs ligand neg
    prot_pos_res = [str(num) + str(res) for res, num in zip(prot_pos.resnames, prot_pos.resnums)]
    for position in lig_neg.positions:
        scores = switch(np.linalg.norm(prot_pos.positions - position, axis=1))
    
    salt_data = diction_me(prot_pos_res, scores, salt_data)

    # protein neg vs ligand pos
    prot_neg_res = [str(num) + str(res) for res, num in zip(prot_neg.resnames, prot_neg.resnums)]
    for position in lig_pos.positions:
        scores = switch(np.linalg.norm(prot_neg.positions - position, axis=1))
    
    salt_data = diction_me(prot_neg_res, scores, salt_data)

    return salt_data

def calculate_polar_interaction(u, ligand, polar_data = {}):
    """
    Calculate polar interaction scores between a ligand and surrounding protein atoms.

    Args:
        u (MDAnalysis.Universe): MDAnalysis Universe object representing the molecular system.
        ligand (str): Selection string for the ligand of interest.
        polar_data (dict, optional): Dictionary to store polar interaction scores (default is an empty dictionary).

    Returns:
        dict: Updated dictionary of polar interaction scores, where each protein residue is associated with its cumulative score.
    """

    lig_back_polar = u.select_atoms(f"{ligand} and (backbone) and (name O* N*)")
    lig_polar = u.select_atoms(f"{ligand} and (resname THR SER GLN ASN CYS MET HSD TRP TYR) and (name O* N* S*)")
    lig_polar += lig_back_polar
    prot_pos = u.select_atoms(f"(around 7 {ligand}) and (resname ARG LYS) and (name CZ NZ)")
    
    # protein pos vs ligand neg
    prot_pos_res = [str(num) + str(res) for res, num in zip(prot_pos.resnames, prot_pos.resnums)]
    
    for position in lig_polar.positions:
        scores = switch(np.linalg.norm(prot_pos.positions - position, axis=1))
    polar_data = diction_me(prot_pos_res, scores, polar_data)

    # protein neg vs ligand pos
    lig_pos = u.select_atoms(f"{ligand} and (resname ARG LYS) and (name CZ NZ)")
    prot_back = u.select_atoms(f"(around 7 {ligand}) and (backbone) and (name O* N*)")
    prot_polar = u.select_atoms(f"(around 7 {ligand}) and (resname THR SER GLN ASN CYS MET HSD TRP TYR) and (name O* N* S*)")
    prot_polar += prot_back
    prot_polar_res = [str(num) + str(res) for res, num in zip(prot_polar.resnames, prot_polar.resnums)]
    
    for position in lig_pos.positions:
        scores = switch(np.linalg.norm(prot_polar.positions - position, axis = 1))
    polar_data = diction_me(prot_polar_res, scores, polar_data)

    return polar_data

def calculate_h_bond(u, ligand, step, HBond = {}):
    """
    Calculate hydrogen bond interactions between a ligand and surrounding atoms.

    Args:
        u (MDAnalysis.Universe): MDAnalysis Universe object representing the molecular system.
        ligand (str): Selection string for the ligand of interest.
        step (int): Number of frames to skip.
        HBond (dict, optional): Dictionary to store hydrogen bond interaction scores (default is an empty dictionary).

    Returns:
        dict: Updated dictionary of hydrogen bond interaction scores, where each residue involved in a hydrogen bond is associated with its cumulative score.
    """

    hbondsa = HBA(universe=u)
    hbondsa.hydrogens_sel = f"(around 7 {ligand}) and (name H* and bonded name O* N*)"
    hbondsa.acceptors_sel = f"{ligand} and (name O* N*)"
    hbondsa.run(step=step)
    resulta = hbondsa.count_by_ids()
    hb_a = [str(num) + str(res) for res, num in zip(u.atoms[resulta[:,1]].resnames, u.atoms[resulta[:,1]].resnums)]

    hbondsb = HBA(universe=u)
    hbondsb.hydrogens_sel = f"{ligand} and (name H* and bonded name O* N*)"
    hbondsb.acceptors_sel = f"(around 7 {ligand}) and (name O* N*)"
    hbondsb.run(step=step)
    resultb = hbondsb.count_by_ids()
    hb_b = [str(num) + str(res) for res, num in zip(u.atoms[resulta[:,2]].resnames, u.atoms[resulta[:,2]].resnums)]

    HBond = diction_me(hb_a, resulta[:,3], HBond)
    HBond = diction_me(hb_b, resultb[:,3], HBond)
    
    return HBond

def calculate_p_stacking(u, ligand, PStack = {}):
    """
    Calculate pi-stacking interactions between a ligand and surrounding aromatic residues in the protein.

    Args:
        u (MDAnalysis.Universe): MDAnalysis Universe object representing the molecular system.
        ligand (str): Selection string for the ligand of interest.
        PStack (dict, optional): Dictionary to store pi-stacking interaction scores (default is an empty dictionary).

    Returns:
        dict: Updated dictionary of pi-stacking interaction scores, where each aromatic residue involved in pi-stacking is associated with its cumulative score.
    """
    
    tyrphe_lig = u.select_atoms(f"{ligand} and resname TYR PHE and name CG CE1 CE2 ")
    his_lig = u.select_atoms(f"{ligand} and resname HIS and name CG CE1 NE2 ")
    trp_lig = u.select_atoms(f"{ligand} and resname TRP and name CG CE2 CE3 ")
    arom_lig = tyrphe_lig + his_lig + trp_lig

    tyrphe_prot = u.select_atoms(f"(resname TYR PHE) and (name CG CE1 CE2)")
    his_prot = u.select_atoms(f"(resname HIS) and (name CG CE1 NE2)")
    trp_prot = u.select_atoms(f"(resname TRP) and (name CG CE2 CE3)")
    arom_prot = tyrphe_prot + his_prot + trp_prot - arom_lig

    arom_lig3 = [arom_lig[i:i + 3] for i in range(0, len(arom_lig), 3)]
    arom_prot3 = [arom_prot[i:i + 3] for i in range(0, len(arom_prot), 3)]

    for lig_arom in arom_lig3:
        for prot_arom in arom_prot3:
            c1 = calculate_center(lig_arom.positions)
            c2 = calculate_center(prot_arom.positions)
            distance = np.linalg.norm(c1 - c2)
            if distance > 6:
                continue

            # Calcola l'angolo di incidenza tra i due piani
            angle = calculate_angle_between_planes(lig_arom.positions, prot_arom.positions)
            if angle > 75 and angle < 105:
                res = str(prot_arom.resnum) + str(prot_arom.resname)
                if res in PStack:
                    PStack[res] += 1
                else:
                    PStack[res] = 1
            elif angle > 165 and angle < 195:
                res = str(prot_arom.resnum) + str(prot_arom.resname)
                if res in PStack:
                    PStack[res] += 1
                else:
                    PStack[res] = 1

    return PStack

def protein_contacts(args: argparse.ArgumentParser):
    """
    Calculate protein-protein interactions and generate visualizations.

    Args:
        args (argparse.ArgumentParser): Command-line arguments containing:
            select (str): Selection string for the ligand of interest.
            topol (file): Topology file for reference (.tpr, .pdb, .psf).
            trj (file): Trajectory file (.trr, .xtc, .dcd).
            step (int): Number of frames to skip.
            limit (float): Threshold value for filtering interactions.

    Returns:
        None: Saves interaction data and visualizations to files.
    """

    coord_dfs = []
    ligand = args.select
    step = args.step

    u = MDAnalysis.Universe(args.topol, args.trj)
    for ts in u.trajectory[:step]:

        hydrophobic = calculate_hydrophobic_interaction(u, ligand)
        salt_bridge = calculate_salt_bridge_interaction(u, ligand)
        polar = calculate_polar_interaction(u, ligand)
        pstack = calculate_p_stacking(u, ligand)
    hbond = calculate_h_bond(u, ligand, step)

    for interaction in [hydrophobic, salt_bridge, polar, pstack]:
        if len(interaction) == 0:
            interaction["0GLY"] = 0 
        df = pd.DataFrame.from_dict(interaction, orient="index")
        df.columns = ["val"]
        df["val"] = df["val"] / df["val"].max()
        coord_dfs.append(df)

    h = pd.DataFrame.from_dict(hbond, orient="index")
    h.columns = ["val"]
    h["val"] = h["val"] / h["val"].max()
    coord_dfs.append(h)

    df_all = pd.concat(coord_dfs, axis = 1)
    df_all = df_all.fillna(0)
    df_all = df_all[(df_all.T > args.limit).any()]
    df_all = df_all.sort_index(ascending = False)
    df_all.columns = ["Hydrophobic", "SaltBridge", "Polar", "p Stack", "H-bond"]
    df_all.to_csv(f"coord_{ligand}.csv")
    ax = df_all.plot.bar(stacked=True, color=["purple", "mediumblue", "lightgreen", "black", "aqua"])
    ax.set_ylim(top=4)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="best")
    ax.set_ylabel("coordination")
    (ax.figure).savefig(f"coord_{ligand}.png", format="png", bbox_inches="tight")


if __name__ == "__main__":
    parser = _parse_CLAs()
    args = parser.parse_args()

    protein_contacts(args)
