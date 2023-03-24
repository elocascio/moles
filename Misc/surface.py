from __future__ import print_function
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import sys, os
import numpy as np
import platform
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy.interpolate import griddata
#import moviepy.editor as mpy
from matplotlib import colors

fes = sys.argv[1]

data = np.loadtxt(fes)

np.savetxt("festrim",data[:,:3])

fes = np.genfromtxt("festrim")
x = fes[:,0]
y = fes[:,1]
z = fes[:,2]

xcli = np.linspace(min(x), max(x), 100)
ycli = np.linspace(min(y), max(y), 100)
X, Y = np.meshgrid(xcli, ycli)
# GRID DATA METHODS: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
Z = griddata((x, y), z, (X, Y), method='linear')

plt.figure(figsize=(11, 10))
plt.imshow(Z, extent=[min(x),max(x),min(y),max(y)], aspect='auto',origin='lower', cmap = 'jet', interpolation='bicubic', alpha=1)
plt.colorbar()
CS = plt.contour(X, Y, Z, 15, colors='k', alpha=0.7, linewidths=0.8, linestyles='dashed', extend='neither')
plt.clabel(CS)
plt.savefig("surface.png", format = "png", dpi = 600)
