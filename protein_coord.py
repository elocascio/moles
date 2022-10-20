def switch(r, r_0=6,a=6,b=12):
    return (1-(r/r_0)**a)/(1-(r/r_0)**b)

def switch_functions(r, r0, r_0 = 6, a = 6, b = 12):
    r = np.asarray(r)
    return np.sum((1-(r/r_0)**a)/(1-(r/r_0)**b))


import MDAnalysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools

def protein_contacts(topol = 'MD.pdb', trj = 'MD.xtc', step = 10):

    Hydrophobic = pd.DataFrame(); Hydrophobic_time = pd.DataFrame()
    PiStacking_T = pd.DataFrame()
    PiStacking_P = pd.DataFrame()
    PiCation = pd.DataFrame()
    SaltBridge = pd.DataFrame(); SaltBridge_time = pd.DataFrame()
    HBond = []
    coord_dfs = []

    u = MDAnalysis.Universe(topol, trj)



    for ts in u.trajectory[::step]:

########### Hydrophobic ###############    
        df_H = pd.DataFrame()
        lig_c  = u.select_atoms(f"segid {ligand} and (resnum 7-14) and (name C*)")
        prot_c = u.select_atoms(f"segid PROA and (around 7 segid {ligand}) and name C*")
        gen = itertools.cycle(range(len(lig_c)))
        prot_list = []
        for res,num in list(zip(prot_c.resnames, prot_c.resnums)):
            prot_list.append(str(num)+str(res))
        for position in lig_c.positions:
            scores = switch(np.linalg.norm(prot_c.positions - position, axis = 1))
            df_H[next(gen)] = scores
        df_H = pd.DataFrame((prot_list, df_H.sum(axis = 1))).T
        df_H[1] = df_H[1].astype(float)
        df_H = df_H.groupby(0).mean()
        Hydrophobic = pd.concat([Hydrophobic, df_H])
        Hydrophobic_time = pd.concat([Hydrophobic_time, df_H], axis = 1)
########### Hydrophobic ###############    
