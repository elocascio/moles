import numpy as np 
import matplotlib.pyplot as plt
import base64
from io import BytesIO

def plot_xvg(xvg, title, filename, xlabel = 'x', ylabel = 'y'):
    a = []
    b = []
    with open(xvg) as xvg:
        for line in xvg.readlines():
            if line.startswith('#') or line.startswith('@'):
                pass
            else:
                x = line.split()
                a.append(float(x[0]))
                b.append(float(x[1]))
        a = np.array(a)
        b = np.array(b)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(a / 1000, b)
        plt.savefig(filename, format = 'png', dpi = 600)        
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        plt.close()
        figfile.seek(0)
        figdata_png = base64.b64encode(figfile.getvalue()).decode()
        plot_string = f'<img src="data:image/png;base64,{figdata_png}" /> '
    return a, b, plot_string

def detachmet(contacts):
    contacts_mean = np.mean(contacts)
    if (np.count_nonzero(contacts < 150) / contacts.shape[0]) < 0.40:
        status = 'poor ligand'
        pass
    else: status = 'attached'
    
    return status, contacts_mean


import yagmail

def send_mail(destination, subject, content, attachment):
    yag = yagmail.SMTP(user = "dinamica.molecolare@gmail.com", password = "AutoNotif1!")
    yag.send(
        to = destination,
        subject = subject,
        contents = content,
        attachments = attachment,
    )
    print('email sent')

import nvsmi
def gpu_manager():
    ids = []
    for process in nvsmi.get_gpu_processes():
        gpu_id = int(str(process).split('gpu_id: ')[-1].split(' | ')[0])
        ids.append(gpu_id)
    print(ids, ids.count(1), ids.count(0))
    if len(ids) > 0:
        if ids.count(0) == 0:
            return 0
        elif ids.count(1) / ids.count(0) > 1 or ids.count(0) == 0:
            return 0
        elif ids.count(1) / ids.count(0) < 1 or ids.count(1) == 0:
            return 1
        else:
            return np.random.randint(2)
    else:
        return np.random.randint(2)

import MDAnalysis
from MDAnalysis.analysis import contacts
from io import BytesIO
import base64

def switch_function(r, r0, r_0 = 6, a = 6, b = 12):
    r = np.asarray(r)
    return np.sum((1-(r/r_0)**a)/(1-(r/r_0)**b))

def coordination():
    u = MDAnalysis.Universe('MD.tpr', 'MD.xtc')
    u_pdb = MDAnalysis.Universe('MD.pdb')
    pdb = []
    for res in u_pdb.segments[0].residues: pdb.append(str(res.resnum) + str(res.resname))
    tpr = []
    for res in u.segments[0].residues: tpr.append(str(res.resnum) + str(res.resname))
    Rosetta = dict(zip(tpr,pdb)) # Rosetta['72PHE'] -> 318PHE
    # Hydrophobic
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
    
    #Polar
    P = []
    p_lig = "(segid seg_1*) and (type O* N* S*)"
    lig_p = u.select_atoms(p_lig)
    prot_p = u.select_atoms("(around 10 segid seg_1*) and (name N* O* S*) and not (resname SOL)")
    for at in prot_p:
        resid = at.resid
        name = at.name
        prot = f"(around 10 segid seg_1*) and (resid {resid}) and (name {name})"
        sel = u.select_atoms(prot)
        coord_P = contacts.Contacts(u, select = (p_lig, prot), refgroup = (lig_p, sel), method=switch_function, kwargs={'r_0':2.5, 'a':8, 'b':12}).run(step=10)
        P.append([str(at.residue.resnum) + str(at.residue.resname), np.mean(coord_P.timeseries[:, 1])])
    
    df_C = pd.DataFrame(H, columns=['residue', 'coord_C'])
    df_P = pd.DataFrame(P, columns=['residue', 'coord_P'])
    aggregation_functions = {'coord_C': 'sum', 'coord_P': 'sum'}
    df_C = df.groupby(df_C['residue']).aggregate(aggregation_functions)
    df_P = df.groupby(df_P['residue']).aggregate(aggregation_functions)
    df_all = pd.concat([df_C, df_P], axis = 1)
    df_all = df_all.fillna(0) ; df_all = df_all[(df_all.T != 0).any()]

    fig, ax = plt.subplots()
    ax.bar(df_all.index, df_all.coord_C)
    ax.bar(df_all.index, df_all.coord_P)

    figfile = BytesIO()
    fig.savefig(figfile, format='png')
    figfile.seek(0)
    figdata_png = base64.b64encode(figfile.getvalue()).decode()
    coord_plot = f'<img src="data:image/png;base64,{figdata_png}" /> '

    return coord_plot
