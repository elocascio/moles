import numpy as np 
import matplotlib.pyplot as plt

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
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.plot(a / 1000, b)
        fig.savefig(filename, format = 'png', dpi = 600)
        plt.close()
    return a, b, fig

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
        if ids.count(1) == 0:
            return 1
        elif ids.count(0) / ids.count(1) > 1 or ids.count(1) == 0:
            return 1
        elif ids.count(0) / ids.count(1) < 1 or ids.count(0) == 0:
            return 0
        else:
            return np.random.randint(2)
    else:
        return np.random.randint(2)

from io import BytesIO
import base64


def fig2html(fig):
    figfile = BytesIO()
    fig.savefig(figfile, format = 'png', bbox_inches = 'tight')
    figfile.seek(0)
    figdata_png = base64.b64decode(figfile.getvalue()).decode()
    imgstr = f'<img src="data:image/png;base64,{figdata_png}" />'
    return imgstr
