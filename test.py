from classLFexploreOld import *
from classTriangulationLF import *
from classArbre import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
xmin = 0.5
xmax = 2
nx = 33
ny = 33

#gridtype = raw_input("Enter grid type ('log' or 'lin'):\n")
#mtype=raw_input("Enter measure type:\n")
#plottype=raw_input("Enter plot type:\n")

gridtype = 'lin'
mtype = '4vn' # measure type
plottype = 'err'


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
    xinit = np.linspace(xmin,xmax,2)
    yinit = np.linspace(ymin,ymax,2)

    ymin2 = ymin + 0.4*(y[1]-y[0])
    ymax2 = ymax + 0.4*(y[1]-y[0])
    xmin2 = xmin + 0.4*(x[1]-x[0])
    xmax2 = xmax + 0.4*(x[1]-x[0])

    x2 = np.linspace(xmin2,xmax2,nx)
    y2 = np.linspace(ymin2,ymax2,ny)
    x2 = np.delete(x2,-1)
    y2 = np.delete(y2,-1)
    matdir = './Matrices2/eeg_leadfield'
    matdir2 = './MatricesTest2/eeg_leadfield'

xx, yy = np.meshgrid(x, y)
xxinit, yyinit = np.meshgrid(xinit, yinit)
xx = xx.flatten()
yy = yy.flatten()
xxinit = xxinit.flatten()
yyinit = yyinit.flatten()

xx2, yy2 = np.meshgrid(x2, y2)
xx2 = xx2.flatten()
yy2 = yy2.flatten()

# LF matrices computimg/export
test = LF_explore([], x, y)
test2 = LF_explore([], x2, y2)
triang = TriangulationLF(xxinit,yyinit,np.array([[0, 1, 2],[3, 2, 1]]))

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
    triang.addCheckLF(np.array(f['linop']), xx[i], yy[i])

for i in range(len(xx2)):
    cond_file_gener(condFile,[xx2[i],yy2[i]])
    savedir = matdir2 + str(i)+'.mat'
    #gain_eeg = LFcomputing(condFile, geomFile, dipoleFile, electrodesFile,savedir)
    f = h5py.File(savedir)
    test2.addLF(np.array(f['linop']))


threshold = 0.01
norm = 'pinv'

arbre = GridArbre(np.array([x[0],y[0]]),np.array([x[-1],y[-1]]))
Err_rectangle = error_calcul(x,y,triang.checkLF,arbre,threshold,norm)
print 'Number of points (rectangle): ' + str(len(unique_points(arbre.get_points())))

plt.figure()
ax = plt.gca()
rec = arbre.get_rectangles()
for r in rec:
    c = r[0]
    ax.add_patch(patches.Rectangle((c[0], c[1]),r[1],r[2],fill=False))
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
plt.xlabel('Brain')
plt.ylabel('Skull')
plt.show()

xxx = np.reshape(xx,(nx,ny))
yyy = np.reshape(yy,(nx,ny))
Err = np.reshape(Err_rectangle,(nx,ny))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xxx, yyy, Err, rstride=1, cstride=1, cmap=cm.coolwarm)
plt.xlabel('Brain')
plt.ylabel('Skull')
plt.show()


LFtest = test.LFinterpol(xx2,yy2)
triang.addLF(triang.checkLF[0])
triang.addLF(triang.checkLF[nx-1])
triang.addLF(triang.checkLF[-nx])
triang.addLF(triang.checkLF[-1])

Err_triangle = triang.make_grid(threshold,'refine',norm)
print 'Number of points (triangle): ' + str(len(triang.x))
plt.figure()
plt.triplot(triang, 'bo-')
plt.xlabel('Brain')
plt.ylabel('Skull')
plt.show()


#xx = np.reshape(xx,(nx,ny))
#yy = np.reshape(yy,(nx,ny))
Err = np.reshape(Err_triangle,(nx,ny))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xxx, yyy, Err, rstride=1, cstride=1, cmap=cm.coolwarm)
plt.xlabel('Brain')
plt.ylabel('Skull')
plt.show()




zz = test.VTmeasure(mtype) # measure of variation
zz2 = np.array([]); # measue of matrices
zztest = np.array([]);
for i in range(len(xx)):
    zz2 = np.append(zz2,np.linalg.norm(test.LF[i]))
for i in range(len(xx2)):
    zztest = np.append(zztest,np.linalg.norm(LFtest[i] - test2.LF[i])/np.linalg.norm(test2.LF[i]))



# Visu

xx = np.reshape(xx,(nx,ny))
yy = np.reshape(yy,(nx,ny))
xx2 = np.reshape(xx2,(nx-1,ny-1))
yy2 = np.reshape(yy2,(nx-1,ny-1))
zz = np.reshape(zz,(nx,ny))
zz2 = np.reshape(zz2,(nx,ny))
zztest = np.reshape(zztest,(nx-1,ny-1))

if plottype != 'mesh':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if plottype == 'var':
        ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
    if plottype == 'norm':
        ax.plot_surface(xx[1:-1,1:-1], yy[1:-1,1:-1], zz2[1:-1,1:-1], rstride=1, cstride=1, cmap=cm.coolwarm)
    if plottype == 'err':
        ax.plot_surface(xx2, yy2, zztest, rstride=1, cstride=1, cmap=cm.coolwarm)
    plt.xlabel('Brain')
    plt.ylabel('Skull')
    plt.show()
else:
    plt.plot(xx, yy,'ro')
    plt.plot(xx2,yy2,'bx')
    plt.show()

