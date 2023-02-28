def switch(r, r_0=6,a=6,b=12):
    return (1-(r/r_0)**a)/(1-(r/r_0)**b)

from inspect import stack
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
import progressbar
from Misc.moles import init
import pickle as pkl
import itertools

def contacts(pdb = 'MD.pdb', xtc = 'MD.xtc', step = 10, ligand = 'UNK'):
    Hydrophobic = pd.DataFrame()
    PiStacking_T = []
    PiStacking_P = []
    PiCation = []
    HBond = []
    SaltBridge = []
    WaterBridge = []
    metal_coord = []
    coord_dfs = []

    print('contact analysis')
    u = MDAnalysis.Universe(pdb,xtc)
    complesso = u.select_atoms(f'(protein) or (resname {ligand})')
    
    if "TIP3" in u.residues.resnames:
        wat = "TIP3"
    elif "SOL" in u.residues.resnames:
        wat = "SOL"
    else:
        print("can't help you bro!")
        exit()
    print(f'water name: {wat}')

    ############# NEG ATOM SELECTION ###############
    ion = None
    neg_atoms = u.select_atoms('(resname {ligand}) and (name O* N* S* Cl* F* Br* I*)')
    for resname in u.residues.resnames:
        if resname in METAL_IONS:
            ion = resname
            print(f'found {ion}')
    if ion != None:
        ion_group = u.select_atoms(f'resname {ion}')
    ###############################################

    water = u.select_atoms(f'(around 10 resname {ligand}) and (resname {wat})') 
    if ion != None:
        complexo = complesso + water + ion_group; comp_ion = complesso + ion_group
    else:
        complexo = complesso + water
##################################### dictionaries #########################################################
    res_list = []
    if ion != None:
        for numres, nameres in list(zip(comp_ion.resnums, comp_ion.resnames)): res_list.append(str(numres)+ str(nameres))
        dic = {key : [] for key in res_list}
    else:
        for numres, nameres in list(zip(complexo.resnums, complexo.resnames)): res_list.append(str(numres)+ str(nameres))
        dic = {key : [] for key in res_list}
    Hydrophobic_dic, PiStacking_dic, PiCation_dic, HBond_dic, SaltBridge_dic = [dic] * 5
###########################################################################################################
    unk = u.select_atoms(f'resname {ligand}'); lig_id = unk.resids[0]; lig_chain = unk.chainIDs[0]

    counter = 1
    with open("ghost", "w"):
        for ts in u.trajectory:
            if ts.frame % step == 0:
                complexo.write('trajj.pdb')
                mol = PDBComplex()
                mol.load_pdb('trajj.pdb')
                UNK = f'{ligand}:{lig_chain}:{lig_id}'
                mol.analyze()
                my_interactions = mol.interaction_sets[UNK]

############### METAL COORDINATION #################
                if ion != None:
                    un = MDAnalysis.Universe('trajj.pdb')
                    neg_atoms = un.select_atoms(f'(resname {ligand}) and (name O* N* S* Cl* F* Br* I*)')
                    ion_group = un.select_atoms(f'resname {ion}')
                    metal_distance = np.linalg.norm(neg_atoms.positions - ion_group.positions, axis = 1)
                    coordination_scores = (metal_distance < 3).astype(int)
                    for coordination_score in coordination_scores:
                        metal_coord.append([str(ion_group.resnames[0]), coordination_score])
####################################################
################# HYDROPHOBIC ########################
                df = pd.DataFrame()
                lig_c  = u.select_atoms(f"resname {ligand} and name C*")
                prot_c = u.select_atoms(f"protein and (around 7 resname {ligand}) and name C*")
                gen = itertools.cycle(range(len(lig_c)))
                prot_list = []
                for res,num in list(zip(prot_c.resnames, prot_c.resnums)):
                    prot_list.append(str(num)+str(res))
                for position in lig_c.positions:
                    scores = switch(np.linalg.norm(prot_c.positions - position, axis = 1))
                    df[next(gen)] = scores
                df = pd.DataFrame((prot_list, df.sum(axis = 1))).T
                df[1] = df[1].astype(float)
                df = df.groupby(0).mean()
                Hydrophobic = pd.concat([Hydrophobic, df])
################# HYDROPHOBIC ########################
#                for hydrophob in my_interactions.hydrophobic_contacts:
#                    Hydrophobic.append([str(hydrophob.resnr) + str(hydrophob.restype), switch(float(hydrophob.distance))])
#                    Hydrophobic_dic[str(hydrophob.resnr) + str(hydrophob.restype)].append(switch(float(hydrophob.distance)))
#                for res in dic.keys():
#                    if len(Hydrophobic_dic[res]) == counter:pass
#                    else: Hydrophobic_dic[res].append(0)
#                    print(Hydrophobic_dic)
############################################################
################ PI STACKING ##############################
                for pi_stacking in my_interactions.pistacking:
                    t = 0; p = 0
                    if pi_stacking.type == 'T':
                        t = switch(float(pi_stacking.distance), r_0=5.5)
                        PiStacking_T.append([str(pi_stacking.resnr) + str(pi_stacking.restype), t])
                    if pi_stacking.type == 'P':
                        p = switch(float(pi_stacking.distance), r_0=5.5)
                        PiStacking_P.append([str(pi_stacking.resnr) + str(pi_stacking.restype), p])
                    
                        stack = p+t
                        PiStacking_dic[str(pi_stacking.resnr) + str(pi_stacking.restype)].append(stack)
                    for res in dic.keys():
                        if len(PiStacking_dic[res]) == counter:pass
                        else: PiStacking_dic[res].append(0)
##########################################################
################ SALT BRIDGE #############################
                for salt_bridge in my_interactions.saltbridge_lneg:
                    salt = switch(float(salt_bridge.distance), r_0=5)
                    SaltBridge.append([str(salt_bridge.resnr) + str(salt_bridge.restype), salt])
                    SaltBridge_dic[str(salt_bridge.resnr) + str(salt_bridge.restype)].append(salt)
                for salt_bridge in my_interactions.saltbridge_pneg:
                    salt = switch(float(salt_bridge.distance), r_0=5)
                    SaltBridge.append([str(salt_bridge.resnr) + str(salt_bridge.restype), salt])
                    SaltBridge_dic[str(salt_bridge.resnr) + str(salt_bridge.restype)].append(salt)
                for res in dic.keys():
                    if len(SaltBridge_dic[res]) == counter:pass
                    else: SaltBridge_dic[res].append(0)
######################################################
################ HBOND ###############################
                for h_bond in my_interactions.hbonds_ldon:
                    HBond.append([str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])
                    if h_bond.type == 'strong': h = 1
                    elif h_bond.type == 'weak': h = 0.5
                    HBond_dic[str(h_bond.resnr) + str(h_bond.restype)].append(h)
                for h_bond in my_interactions.hbonds_pdon:
                    HBond.append([str(h_bond.resnr) + str(h_bond.restype), str(h_bond.type)])
                    if h_bond.type == 'strong': h = 1
                    elif h_bond.type == 'weak': h = 0.5
                    HBond_dic[str(h_bond.resnr) + str(h_bond.restype)].append(h)
                for res in dic.keys():
                    if len(HBond_dic[res]) == counter:pass
                    else: HBond_dic[res].append(0)
#####################################################
################ PI CATION ##########################
                for picat in my_interactions.pication_laro:
                    PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])
                    PiCation_dic[str(picat.resnr) + str(picat.restype)].append(switch(float(picat.distance)))
                for picat in my_interactions.pication_paro:
                    PiCation.append([str(picat.resnr) + str(picat.restype), switch(float(picat.distance), r_0=4)])
                    PiCation_dic[str(picat.resnr) + str(picat.restype)].append(switch(float(picat.distance)))
                for res in dic.keys():
                    if len(PiCation_dic[res]) == counter:pass
                    else: PiCation_dic[res].append(0)
####################################################
################ WATER BRIDGE ######################                
                for water_bridge in my_interactions.water_bridges:
                    WaterBridge.append([str(water_bridge.resnr) + str(water_bridge.restype), water_bridge.type])
####################################################

                counter += 1
    Hydrophobic = Hydrophobic.reset_index()
    Hydrophobic = Hydrophobic.values
    
#    n = 0
#    for dict in [Hydrophobic_dic, PiStacking_dic, PiCation_dic, HBond_dic, SaltBridge_dic]:
#        print(n)
#        dictframe = pd.DataFrame.from_dict(dict); dictframe.to_csv(f"{n}.csv")
    

    for coord, coord_type in list(zip([Hydrophobic, PiStacking_T, PiStacking_P, SaltBridge, HBond, PiCation, WaterBridge, metal_coord], ['Hydrophobic', 'Pi_stacking_T', 'Pi_stacking_P', 'Salt_Bridge', 'H_bond', 'Pi_Cation', 'Water_Bridge', 'Metal_Coordination'])):
        
        if coord in [HBond, WaterBridge]:
            df = pd.DataFrame(coord, columns=['residue', coord_type])
            aggregation = {coord_type: 'count'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            if args.abs:
                df[coord_type] = df[coord_type].apply(lambda x: x / (u.trajectory.n_frames / step))
                coord_dfs.append(df)
            else:
                df[coord_type] = df[coord_type].apply(lambda x: x / df[coord_type].max())
                coord_dfs.append(df)
        elif coord_type in "Hydrophobic":
            df = pd.DataFrame(coord, columns = ['residue', coord_type])
            aggregation = {coord_type: 'sum'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            if args.abs:
                df[coord_type] = df[coord_type].apply(lambda x: x / (u.trajectory.n_frames / step))
                coord_dfs.append(df)
            else:
                df[coord_type] = df[coord_type].apply(lambda x: x / df[coord_type].max())
                coord_dfs.append(df)
        else: 
            df = pd.DataFrame(coord, columns = ['residue', coord_type])
            aggregation = {coord_type: 'sum'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            if args.abs:
                df[coord_type] = df[coord_type].apply(lambda x: x / (u.trajectory.n_frames / step))
                coord_dfs.append(df)
            else:
                df[coord_type] = df[coord_type].apply(lambda x: x / df[coord_type].max())
                coord_dfs.append(df)
    
    df_all = pd.concat(coord_dfs, axis = 1); df_all = df_all.fillna(0); df_all = df_all[(df_all.T > 0.2).any()]
    df_all = df_all.sort_values(by=['residue'])
    if args.abs:
        df_all.to_csv(f'coord_abs_{ligand}.csv')
    else:
        df_all.to_csv(f'coord_{ligand}.csv')
    ax = df_all.plot.bar(stacked = True, color = colors)
    ax.set_ylim(top=4)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='best')
    ax.set_ylabel('coordination')
    if args.abs:
        (ax.figure).savefig(f'coord_{ligand}_abs.png', format = 'png', bbox_inches = 'tight')
    else:
        (ax.figure).savefig(f'coord_{ligand}.png', format = 'png', bbox_inches = 'tight')

    figfile = BytesIO()
    (ax.figure).savefig(figfile, format='png', bbox_inches = 'tight')
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue()).decode()
    string = f'<img src="data:image/png;base64,{figdata_png}" /> '
    plt.close()
    return string

if __name__=='__main__':
#    init()
    parser = argparse.ArgumentParser()
    parser.add_argument("-lig", type=str, help="name of ligand --- default= UNK", default="UNK")
    parser.add_argument("-pdb", type= str, help="file pdb gromacs --- default MD.pdb", default='MD.pdb')
    parser.add_argument("-xtc", type= str, help="file xtc gromacs --- default MD.xtc", default='MD.xtc')
    parser.add_argument("-step", type= int, help="value of step --- default 10", default='10')
    parser.add_argument("-abs", action= "store_true", help="yield absolute contant values")
    args= parser.parse_args()
    contacts(ligand = args.lig, pdb = args.pdb, xtc = args.xtc, step = args.step)
