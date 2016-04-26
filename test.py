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

"""
ymin = 0
ymax = 2
xmin = 0
xmax = 3
nx = 3
ny = 3
"""
#gridtype = raw_input("Enter grid type ('log' or 'lin'):\n")
#mtype=raw_input("Enter measure type:\n")

gridtype = 'lin'
mtype = '4vn' # measure type


if gridtype == 'log':
    x = np.logspace(np.log10(xmin),np.log10(xmax),nx)
    y = np.logspace(np.log10(ymin),np.log10(ymax),ny)
    matdir = './Matriceslog/eeg_leadfield'
    
    x2, y2 = [], []
    for i in range(len(x)-1):
        x2 = np.append(x2,0.6*x[i]+0.4*x[i+1])
    for i in range(len(y)-1):
        y2 = np.append(y2,0.6*y[i]+0.4*y[i+1])
    matdir2 = './MatriceslogTest/eeg_leadfield'

if gridtype == 'lin':
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)

    ymin2 = ymin + 0.4*(y[1]-y[0])
    ymax2 = ymax + 0.4*(y[1]-y[0])
    xmin2 = xmin + 0.4*(x[1]-x[0])
    xmax2 = xmax + 0.4*(x[1]-x[0])

    x2 = np.linspace(xmin2,xmax2,nx)
    y2 = np.linspace(ymin2,ymax2,ny)
    x2 = np.delete(x2,-1)
    y2 = np.delete(y2,-1)
    matdir = './Matrices/eeg_leadfield'
    matdir2 = './MatricesTest/eeg_leadfield'

xx, yy = np.meshgrid(x, y)
xx = xx.flatten()
yy = yy.flatten()

xx2, yy2 = np.meshgrid(x2, y2)
xx2 = xx2.flatten()
yy2 = yy2.flatten()

# LF matrices computimg/export
test = LF_explore([], x, y)
test2 = LF_explore([], x2, y2)

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
for i in range(len(xx2)):
    cond_file_gener(condFile,[xx2[i],yy2[i]])
    savedir = matdir2 + str(i)+'.mat'
    #gain_eeg = LFcomputing(condFile, geomFile, dipoleFile, electrodesFile,savedir)
    f = h5py.File(savedir)
    test2.addLF(np.array(f['linop']))
LFtest = test.LFinterpol(xx2,yy2)

zz = test.VTmeasure(mtype) # measure of variation

zz2 = np.array([]); # measue of matrices
zztest = np.array([]);
for i in range(len(xx)):
    zz2 = np.append(zz2,np.linalg.norm(test.LF[i]))
for i in range(len(xx2)):
    zztest = np.append(zztest,np.linalg.norm(LFtest[i] - test2.LF[i])/np.linalg.norm(LFtest[i]))



# Visu
xx = np.reshape(xx,(nx,ny))
yy = np.reshape(yy,(nx,ny))
xx2 = np.reshape(xx2,(nx-1,ny-1))
yy2 = np.reshape(yy2,(nx-1,ny-1))
zz = np.reshape(zz,(nx,ny))
zz2 = np.reshape(zz2,(nx,ny))
zztest = np.reshape(zztest,(nx-1,ny-1))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
#ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz2[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
ax.plot_surface(xx2, yy2, zztest, rstride=1, cstride=1, cmap=cm.coolwarm)
plt.xlabel('Brain')
plt.ylabel('Scull')
plt.show()

"""
plt.plot(xx, yy,'ro')
plt.plot(xx2,yy2,'bx')
plt.show()
"""
