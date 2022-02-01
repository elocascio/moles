import numpy as np 
from os import system as sys
import matplotlib.pyplot as plt
import argparse
from utils import clean_xvg

parser = argparse.ArgumentParser()
parser.add_argument("-trr",type=str,help="file trr gromacs --- default MD.trr",default='MD.trr')
parser.add_argument("-tpr",type=str,help="file tpr gromacs --- default MD.tpr",default= 'MD.tpr')
parser.add_argument("-ndx",type=str,help="file index gromacs --- default index.ndx",default='index.ndx')
parser.add_argument("-method",type=str,help="Method for cluster determination: linkage, jarvis-patrick, monte-carlo, diagonalization, gromos --- default gromos",default='gromos')
parser.add_argument("-cutoff",type=float,help="n of cutoff --- default 0.1",default=0.1,)
parser.add_argument("-fit", type=str, help="fit-least square group")
parser.add_argument("-ogroup", type=str,help="ogroup-group output") 

args = parser.parse_args()

sys (f'''gmx cluster -f {args.trr} -s {args.tpr} -n {args.ndx} -sz -cl -method {args.method} -cutoff {args.cutoff} << EOF
{args.fit}
{args.ogroup}
EOF''')

n_clust,freq = clean_xvg('clust-size.xvg')

fig = plt.figure(figsize=(5,7))
plt.title('CLUSTER')
plt.pie(freq)
labels= np.round(freq/freq.sum()*100, 2)
plt.legend(labels,title="CLUSTER %",bbox_to_anchor=(1, .7))
plt.tight_layout()
fig.savefig(f'clust_{args.cutoff}_{args.method}.png',format='png', dpi=1600)