import numpy as np 
import matplotlib.pyplot as plt
import smtplib

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
        plt.close()
    return a, b

def detachmet(contacts):
    contact0 = contacts[0]
    contacts_mean = np.mean(contacts)
    for contact in contacts:
        if contact / contact0 < 0.40:
            status = 'detachment'
            break
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
import numpy as np
def gpu_manager():
    NGPUS = 2
    ids = []
    for process in nvsmi.get_gpu_processes():
        gpu_id = int(str(process).split('gpu_id: ')[-1].split(' | ')[0])
        ids.append(gpu_id)
    print(ids, ids.count(1), ids.count(0))
    if len(ids) > 0:
        if ids.count(0) / ids.count(1) > 1 or ids.count(1) == 0:
            return 1
        elif ids.count(0) / ids.count(1) < 1 or ids.count(0) == 0:
            return 0
        else:
            return np.random.randint(2)