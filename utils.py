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
    contacts_mean = np.mean(contacts)
    for contact in contacts:
        if contact / contact0 < 0.40:
            status = 'detachmet'
            break
        else: status = 'attached'
    
    return status, contacts_mean

def send_mail(destination, content):
    mail = smtplib.SMTP('smtp.gmail.com', 587)
    mail.ehlo()
    mail.starttls()
    mail.login('dinamica.molecolare@gmail.com', 'AutoNotif1!')
    mail.sendmail('dinamica.molecolare@gmail.com', destination, content)
    mail.close()
    print('send')

import sys
import os
import urllib.request
import urllib.parse
import urllib


def get_zinc_smile(zinc_id):
    """
    Gets the corresponding SMILE string for a ZINC ID query from
    the ZINC online database. Requires an internet connection.
    Keyword arguments:
        zinc_id (str): A valid ZINC ID, e.g. 'ZINC00029323'
        backend (str): zinc12 or zinc15
    Returns the SMILE string for the corresponding ZINC ID.
        E.g., 'COc1cccc(c1)NC(=O)c2cccnc2'
    """
    stripped_id = zinc_id.strip('ZINC')

    min_len = 12
    base_path = 'http://zinc15.docking.org/substances/'
    line_lookup = 'id="substance-smiles-field" readonly value="'
    first_linesplit = line_lookup
    second_linesplit = '">'

    while len(stripped_id) < min_len:
        stripped_id = '0' + stripped_id

    smile_str = None
      
    try:
        response = urllib.request.urlopen('{}{}'
                                      .format(base_path, stripped_id))
    except urllib.error.HTTPError:
        print('Invalid ZINC ID {}'.format(zinc_id))
        response = []
    for line in response:
        line = line.decode(encoding='UTF-8').strip()
        if line_lookup in line:
            line = (line.split(first_linesplit)[-1]
                    .split(second_linesplit)[0])
            smile_str = urllib.unquote(line)
            break
    return smile_str

