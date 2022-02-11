def switch(r, r_0=6,a=6,b=12):
    return (1-(r/r_0)**a)/(1-(r/r_0)**b)

import MDAnalysis
from plip.structure.preparation import PDBComplex
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from var import colors, METAL_IONS
import base64
from io import BytesIO
import argparse
from glob import glob

def contacts(pdb = 'MD.pdb', xtc = 'MD.xtc', step = 10, ligand = 'UNK'):
    Hydrophobic = []
    PiStacking_T = []
    PiStacking_P = []
    PiCation = []
    HBond = []
    SaltBridge = []
    WaterBridge = []
    
    coord_dfs = []

    print('contact analysis')
    u = MDAnalysis.Universe(pdb,xtc)
    complexo = u.select_atoms(f'(protein) or (resname {args.lig})')
    
    if "TIP3" in u.residues.resnames:
        wat = "TIP3"
    elif "SOL" in u.residues.resnames:
        wat = "SOL"
    else:
        print("can't help you bro!")
        exit()
    
    ############# NEG ATOM SELECTION ###############
    tpr = glob("*.tpr")
    utpr = MDAnalysis.Universe(tpr[0])
    neg_atoms = utpr.select_atoms('(resname {args.lig}) and (prop charge  < -0.8)')
    for resname in utpr.residues.resnames:
        if resname in METAL_IONS:
            ion = resname

    ###############################################

    water = u.select_atoms(f'(around 10 resname {args.lig}) and (resname {wat})') 
    complexo = complexo + water
    unk = u.select_atoms(f'resname {ligand}'); unk = unk.resids[0]
    for ts in u.trajectory:
        if ts.time % step == 0:
            complexo.write('trajj.pdb')
            mol = PDBComplex()
            mol.load_pdb('trajj.pdb')
            UNK = f'{ligand}:B:{unk}'
            mol.analyze()
            my_interactions = mol.interaction_sets[UNK]

            for hydrophob in my_interactions.hydrophobic_contacts:
                Hydrophobic.append([str(hydrophob.resnr) + str(hydrophob.restype), switch(float(hydrophob.distance))])


            for pi_stacking in my_interactions.pistacking:
                if pi_stacking.type == 'T':
                    PiStacking_T.append([str(pi_stacking.resnr) + str(pi_stacking.restype), switch(float(pi_stacking.distance), r_0=5.5)])
                if pi_stacking.type == 'P':
                    PiStacking_P.append([str(pi_stacking.resnr) + str(pi_stacking.restype), switch(float(pi_stacking.distance), r_0=5.5)])

            for salt_bridge in my_interactions.saltbridge_lneg:
                SaltBridge.append([str(salt_bridge.resnr) + str(salt_bridge.restype), switch(float(salt_bridge.distance), r_0=5)])
            for salt_bridge in my_interactions.saltbridge_pneg:
                SaltBridge.append([str(salt_bridge.resnr) + str(salt_bridge.restype), switch(float(salt_bridge.distance), r_0=5)])

            for h_bond in my_interactions.hbonds_ldon:
                HBond.append([str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])
            for h_bond in my_interactions.hbonds_pdon:
                HBond.append([str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])

            for picat in my_interactions.pication_laro:
                PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])
            for picat in my_interactions.pication_paro:
                PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])

            for water_bridge in my_interactions.water_bridges:
                WaterBridge.append([str(water_bridge.resnr) + str(water_bridge.restype), water_bridge.type])
        
            np.linalg.norm(neg_atoms.positions - ion.position, axis = 1)

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
    
    df_all = pd.concat(coord_dfs, axis = 1); df_all = df_all.fillna(0); df_all = df_all[(df_all.T > 0).any()]
    df_all = df_all.sort_values(by=['residue'])
    df_all.to_csv('coord.csv')
    ax = df_all.plot.bar(stacked = True, color = colors)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='best')
    ax.set_ylabel('coordination')
    (ax.figure).savefig('coord.png', format = 'png', bbox_inches = 'tight')
    
    figfile = BytesIO()
    (ax.figure).savefig(figfile, format='png', bbox_inches = 'tight')
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue()).decode()
    string = f'<img src="data:image/png;base64,{figdata_png}" /> '
    plt.close()
    return string

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-lig", type=str, help="name of ligand --- default= UNK", default="UNK")
    parser.add_argument("-pdb", type= str, help="file pdb gromacs --- default MD.pdb", default='MD.pdb')
    parser.add_argument("-xtc", type= str, help="file xtc gromacs --- default MD.xtc", default='MD.xtc')
    parser.add_argument("-step", type= int, help="value of step --- default 10", default='10')
    args= parser.parse_args()
    contacts(ligand = args.lig, pdb = args.pdb, xtc = args.xtc, step = args.step)
