from utils import send_mail

content = """
This is an auto-generated email. Do not respond to this email address. 
For further information contact ettore.locascio@unicatt.it.

This is HTMD (High Throughput Molecular Dynamics).
Status : running
Ligand processed : 10/1000
Detached ligand : 8

GoodBye"""

send_mail('stefano.dellalonga@univaq.it', content)