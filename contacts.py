import MDAnalysis
from MDAnalysis.analysis import contacts
import numpy as np
import pandas as pd
from rdkit import Chem
from biopandas.mol2 import PandasMol2
from io import BytesIO
import base64
from utils import fig2html


def switch_function(r, r0, r_0 = 6, a = 6, b = 12):
    r = np.asarray(r)
    return np.sum((1-(r/r_0)**a)/(1-(r/r_0)**b))

def contacts(tpr = 'MD.tpr', xtc = 'MD.xtc', pdb = 'MD.pdb', mol2) 
    u = MDAnalysis.Universe(tpr, xtc)
    u_pdb = MDAnalysis.Universe(pdb)
    if mol = Chem.MolFromMol2File(mol2, removeHs = False):
    else:
        return 'ERROR'

    pdb = []
    for res in u_pdb.segments[0].residues: pdb.append(str(res.resnum) + str(res.resname))
    tpr = []
    for res in u.segments[0].residues: tpr.append(str(res.resnum) + str(res.resname))
    Rosetta = dict(zip(tpr,pdb)) # Rosetta['72PHE'] -> 318PHE
#------------------------------- Hydrophobic ----------------------------------
    H = []
    c_lig = "(segid seg_1*) and (type C*)"
    lig_c = u.select_atoms(c_lig)
    prot_c = u.select_atoms("(around 10 segid seg_1*) and (name C*)")
    for at in prot_c:
        resid = at.resid
        name = at.name
        prot = f"(around 10 segid seg_1*) and (resid {resid}) and (name {name})"
        sel = u.select_atoms(prot)
        coord_C = contacts.Contacts(u, select = (c_lig, prot), refgroup = (lig_c, sel), method=switch_function, kwargs={'r_0':6, 'a':6, 'b':12}).run(step=10)
        H.append([str(at.residue.resnum) + str(at.residue.resname) , np.mean(coord_C.timeseries[:, 1])])
    print('H fatto')    
#------------------------------- Polar -------------------------------------------
    P = []
    p_lig = "(segid seg_1*) and (type O* N* S*)"
    lig_p = u.select_atoms(p_lig)
    prot_p = u.select_atoms("(around 10 segid seg_1*) and (name N* O* S*) and not (resname SOL)")
    for at in prot_p:
        resid = at.resid
        name = at.name
        prot = f"(around 10 segid seg_1*) and (resid {resid}) and (name {name})"
        sel = u.select_atoms(prot)
        coord_P = contacts.Contacts(u, select = (p_lig, prot), refgroup = (lig_p, sel), method=switch_function, kwargs={'r_0':2.5, 'a':8, 'b':12}).run(step=100)
        P.append([str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
    print('P fatto')

#-------------------------------- H bonds ---------------------------------------------
    hbonds = HBA(universe=u)
    unk = "resname UNK"
    protein = "(around 10 resname UNK) and (protein)"
    # ligand hydrogen
    hbonds_list = []
    if len(hbonds.guess_hydrogens(unk)) != 0:
        hbonds.hydrogens_sel = hbonds.guess_hydrogens(unk)
        hbonds.acceptors_sel = hbonds.guess_acceptors(protein)
        unk_h = hbonds.run(step = 100)
        for atoms, count in list(zip(unk_acc.count_by_ids()[:, 2], unk_acc.count_by_ids()[:, 3])):
            at = u.select_atoms(f'index {atoms}')
            if at.resnames[0] == 'SOL':
                pass
            else:
                hbonds_list.append(str(at.resnums[0]) + str(at.resnames[0]), count)
                print([str(at.resnums[0]) + str(at.resnames[0]), count])
    else: pass
    if len(hbonds.guess_acceptors(unk)) != 0:
        hbonds.hydrogens_sel = hbonds.guess_hydrogens(protein)
        hbonds.acceptors_sel = hbonds.guess_acceptors(unk)
        unk_acc = hbonds.run(step = 100)
        for atoms, count in list(zip(unk_acc.count_by_ids()[:, 0], unk_acc.count_by_ids()[:, 3])):
            at = u.select_atoms(f'index {atoms}')
            if at.resnames[0] == 'SOL':
                pass
            else:
                hbonds_list.append([str(at.resnums[0]) + str(at.resnames[0]), count])
                print(str(at.resnums[0]) + str(at.resnames[0]), count)
    else: pass

#----------------------------- Salt Bridge, Pi-stacking and Pi-Cation ----------------------------------
    pmol = PandasMol2().read_mol2(mol2)
    charged_atoms = []
    aromatic_atoms = []

    for x in mol.GetAtoms():
        if x.GetFormalCharge() != 0:
            charged_atoms.append(pmol.df[pmol.df['atom_id'] == x.GetIdx() + 1].atom_name.values[0])
        else: pass
        if x.GetIsAromatic() == True:
            aromatic_atoms.append(pmol.df[pmol.df['atom_id'] == x.GetIdx() + 1].atom_name.values[0])
        else: pass

    # Salt Bridge
    if len(charged_atoms) > 0:
        S = []
        s_lig = f'(resname UNK) and (name {" ".join(charged_atoms)})'
        lig_s = u.select_atoms(s_lig)
        prot_s = u.select_atoms(f'(around 10 resname UNK) and (name NZ NH* OD*)')
        for at in prot_s:
            resid = at.resid
            name = at.name
            prot = f'(around 10 resname UNK) and (resid {resid}) and (name {name})'
            sel = u.select_atoms(prot)
            Coord_S = contacts.Contacts(u, select = (s_lig, prot), refgroup = (lig_s, sel), method = switch_function, kwargs = {'r_0': 5, 'a': 6, 'b': 12}).run(step = 100)
            S.append([str(at.residue.resnum) + str(at.residue.resname), np.mean(Coord_S.timeseries[:, 1])])
        print('S fatto')
    else : pass

    # Aromatic
    if len(aromatic_atoms) >0:
        PS = []
        ps_lig = f'(resname UNK) and (name {" ".join(aromatic_atoms)})'
        lig_ps = u.select_atoms(ps_lig)
        prot_ps = u.select_atoms('(around 10 resname UNK) and (resname PHE TYR TRP) and not (name N H* C O* CA)')
        for at in prot_ps:
            resid = at.resid
            name = at.name
            prot = f'(around 10 resname UNK) and (resid {resid}) and (name {name})'
            sel = u.select_atoms(prot)
            Coord_PS = contacts.Contacts(u, select = (ps_lig, prot), refgroup = (lig_ps, sel), method = switch_function, kwargs = {'r_0': 5.5, 'a': 6, 'b': 12}).run(step = 100)
            PS.append([str(at.residue.resnum) + str(at.residue.resname), np.mean(Coord_PS.timeseries[:, 1])])
        print('PS fatto')
    # Pi-Cation
        PC = []
        prot_pc = u.select_atoms((f'(around 10 resname UNK) and (name NZ NH*)'))
        for at in prot_pc:
            resid = at.resid
            name = at.name
            prot = f'(around 10 resname UNK) and (resid {resid}) and (name {name})'; sel = u.select_atoms(prot)
            Coord_PC = contacts.Contacts(u, select = (ps_lig, prot), refgroup = (lig_ps, sel), method = switch_function, kwargs = {'r_0': 4, 'a': 6, 'b': 12}).run(step = 100)
            PC.append([str(at.residue.resnum) + str(at.residue.resname), np.mean(Coord_PC.timeseries[:, 1])])
            print('PC fatto')
    else: pass

    coord_dfs = []
    for coord,coord_type in list(zip([H, P, hbonds_list, S, PS, PC], ['Hydrophobic', 'Polar', 'H_bonds', 'Salt_Bridge', 'Pi_Stacking', 'Pi_Cation'])):
        try:
            df = pd.DataFrame(coord, columns=['residue', coord_type]); df.residue = df['residue'].replace(Rosetta, regex = True)
            aggregation = {coord_type: 'sum'}
            df = df.groupby(df['residue']).aggregate(aggregation)
            coord_dfs.append(df)
        except:
            print('Empy Dataframe: ', coord, coord_type)
    df_all = pd.concat(coord_dfs, axis = 1); df_all = df_all.fillna(0); df_all = df_all[(df_all.T > 0).any()]
    fig = df_all.plot.bar(stacked = True).figure
    fig.savefig('contacts.png', format = 'png', bbox_inches = 'tight')
    
    return fig
    
    