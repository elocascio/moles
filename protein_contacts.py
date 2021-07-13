def switch(r, r_0=6,a=6,b=12):
    return (1-(r/r_0)**a)/(1-(r/r_0)**b)

def switch_functions(r, r0, r_0 = 6, a = 6, b = 12):
    r = np.asarray(r)
    return np.sum((1-(r/r_0)**a)/(1-(r/r_0)**b))

import MDAnalysis
from MDAnalysis.analysis import contacts
from plip.structure.preparation import PDBComplex
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import base64
from io import BytesIO

def protein_contacts(pdb = 'MD.pdb', xtc = 'MD.xtc', step = 10):
    Hydrophobic = []
    PiStacking_T = []
    PiStacking_P = []
    PiCation = []
    HBond = []
    SaltBridge = []
    WaterBridge = []
    Polar = []
    coord_dfs = []

    u = MDAnalysis.Universe(pdb, xtc)
    complexo = u.select_atoms('(protein) or (resname SOL)')
    
    for ts in u.trajectory:
        if ts.time % step == 0:
            complexo.write('trajj.pdb')
            mol = PDBComplex(); mol.load_pdb('trajj.pdb')
            mol.analyze()
            my_interactions = mol.interaction_sets[list(mol.interaction_sets.keys())[0]]

            for hydrophob in my_interactions.hydrophobic_contacts:
                Hydrophobic.append([str(hydrophob.reschain) + str(hydrophob.resnr) + str(hydrophob.restype), switch(float(hydrophob.distance))])


            for pi_stacking in my_interactions.pistacking:
                if pi_stacking.type == 'T':
                    PiStacking_T.append([str(pi_stacking.reschain) + str(pi_stacking.resnr) + str(pi_stacking.restype), switch(float(pi_stacking.distance), r_0=5.5)])
                if pi_stacking.type == 'P':
                    PiStacking_P.append([str(pi_stacking.reschain) + str(pi_stacking.resnr) + str(pi_stacking.restype), switch(float(pi_stacking.distance), r_0=5.5)])

            for salt_bridge in my_interactions.saltbridge_lneg:
                SaltBridge.append([str(salt_bridge.reschain) + str(salt_bridge.resnr) + str(salt_bridge.restype), switch(float(salt_bridge.distance), r_0=5)])
            for salt_bridge in my_interactions.saltbridge_pneg:
                SaltBridge.append([str(salt_bridge.reschain) + str(salt_bridge.resnr) + str(salt_bridge.restype), switch(float(salt_bridge.distance), r_0=5)])

            for h_bond in my_interactions.hbonds_ldon:
                HBond.append([str(h_bond.reschain) + str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])
            for h_bond in my_interactions.hbonds_pdon:
                HBond.append([str(h_bond.reschain) + str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])

            for picat in my_interactions.pication_laro:
                PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])
            for picat in my_interactions.pication_paro:
                PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])
   
    for coord, coord_type in list(zip([Hydrophobic, PiStacking_T, PiStacking_P, SaltBridge, HBond, PiCation, WaterBridge], ['Hydrophobic', 'Pi_stacking_T', 'Pi_stacking_P', 'Salt_Bridge', 'H_bond', 'Pi_Cation', 'Water_Bridge'])):
        if coord in [HBond, WaterBridge]:
            df = pd.DataFrame(coord, columns=['residue', coord_type])
            aggregation = {coord_type: 'count'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            print(df)
            df[coord_type] = df[coord_type].apply(lambda x: x / (len(u.trajectory) / step))
            coord_dfs.append(df)
        else: 
            df = pd.DataFrame(coord, columns = ['residue', coord_type])
            aggregation = {coord_type: 'sum'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            df[coord_type] = df[coord_type].apply(lambda x: x / (len(u.trajectory) / step))
            coord_dfs.append(df)
        
    # POLAR CONTACTS
    pep_polar_str ='(segid D) and not (backbone) and (name O* N* S*) and not (resname ARG LYS ASP GLU)' ; pep_polar = u.select_atoms(pep_polar_str)
    pep_cation_str = '(segid D) and not (backbone) and (name N*) and (resname ARG LYS )'; pep_cation = u.select_atoms(pep_cation_str)
    protein_polar_str = '(around 10 segid D) and (name O* N* S*) and (resname THR SER ASN GLN CYS MET)'; protein_polar = u.select_atoms(protein_polar_str)
    protein_cation_str = '(around 10 segid D) and not (backbone) and (name N*) and (resname ARG LYS)' ; protein_cation = u.select_atoms(protein_cation_str) 

    P = []
    for at in protein_polar:
        prot = f"(resid {at.resid}) and (name {at.name} and (segid {at.chainID}))"
        coord_P = contacts.Contacts(u,
                                     select = (pep, prot) ,
                                     refgroup = (pep_plus,u.select_atoms(prot)) ,
                                     radius = 100 ,
                                     method=switch_functions ,
                                     kwargs={'r_0':2.5, 'a':8, 'b':12}
                                     ).run(step=100)
    P.append([str(at.chainID) + str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
