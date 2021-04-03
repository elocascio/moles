import numpy as np 
import matplotlib.pyplot as plt
import smtplib

def plot_xvg(xvg, title, xlabel = 'x', ylabel = 'y'):
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
        plt.savefig('contacts.png', format = 'png', dpi = 600)
        plt.close()
    return a, b

def detachmet(contacts):
    contact0 = contacts[0]
    for contact in contacts:
        if contact / contact0 < 0.40:
            status = 'detachmet'
            break
        elif contact / contact0 > 1:
            status = 'incredible ligand, or something suspicious'
        else: status = 'attached'
    
    return status

def send_mail(destination, content):
    mail = smtplib.SMTP('smtp.gmail.com', 587)
    mail.ehlo()
    mail.starttls()
    mail.login('dinamica.molecolare@gmail.com', 'AutoNotif1!')
    mail.sendmail('dinamica.molecolare@gmail.com', destination, content)
    mail.close()
    print('send')