from classLFexplore import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import random
import h5py

import openmeeg as om

# Grid initialization
ymin = 0.01
ymax = 0.1
xmin = 0.1
xmax = 1
nx = 25
ny = 25

gridtype = 'log'
mtype = '4vn' # measure type

if gridtype is 'log':
    x = np.logspace(np.log10(xmin),np.log10(xmax),nx)
    y = np.logspace(np.log10(ymin),np.log10(ymax),ny)
    matdir = './Matriceslog/eeg_leadfield'
if gridtype is 'lin':
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)
    matdir = './Matrices/eeg_leadfield'

xx, yy = np.meshgrid(x, y)
xx = xx.flatten()
yy = yy.flatten()

# LF matrices computimg/export
test = LF_explore([], x, y)

geomFile                    = './data/Head1.geom'
condFile                    = './data/Head1test.cond'
dipoleFile                  = './data/Head1.dip'
electrodesFile              = './data/Head1new.txt'

for i in range(len(xx)):
    cond_file_gener(condFile,[xx[i],yy[i]])
    savedir = matdir + str(i)+'.mat'
    #gain_eeg = LFcomputing(condFile, geomFile, dipoleFile, electrodesFile,savedir)
    f = h5py.File(savedir)
    test.addLF(np.array(f['linop']))

zz = test.VTmeasure(mtype) # measure of variation

zz2 = np.array([]); # measue of matrices
for i in range(len(xx)):
    zz2 = np.append(zz2,np.linalg.norm(test.LF[i]))


# Visu
xx = np.reshape(xx,(nx,ny))
yy = np.reshape(yy,(nx,ny))
zz = np.reshape(zz,(nx,ny))
zz2 = np.reshape(zz2,(nx,ny))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz2[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
plt.xlabel('Brain')
plt.ylabel('Scull')
plt.show()


