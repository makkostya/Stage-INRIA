import numpy as np
import math
import sys
import openmeeg as om

def modulo(i, imax):
    return min(abs(i),imax) - max(0,i-imax)
def cond_file_gener(filename, sigma):
    lines = []
    with open(filename) as infile:
	i=0
        for line in infile:
	    if i==4 or i==5:
		r = line.split()
            	line = line.replace(r[1], str(sigma[i-4]))
	    i = i+1
            lines.append(line)
    with open(filename, 'w') as outfile:
        for line in lines:
            outfile.write(line)

def LFcomputing(condFile, geomFile, dipoleFile, electrodesFile,savedir):
    """
    condFile = 'om_demo.cond'
    geomFile = 'om_demo.geom'
    dipoleFile = 'cortex.dip'
    squidsFile = 'meg_squids.txt'
    electrodesFile = 'eeg_electrodes.txt' 
    """
    # Load data
    geom = om.Geometry()
    geom.read(geomFile,condFile)
    dipoles = om.Matrix()
    dipoles.load(dipoleFile)
    #squids = om.Sensors()
    #squids.load(squidsFile)
    electrodes = om.Sensors()
    electrodes.load(electrodesFile)

    # Compute forward problem
    gaussOrder = 3
    use_adaptive_integration = True

    hm = om.HeadMat(geom,gaussOrder)
    hminv = hm.inverse()
    dsm = om.DipSourceMat (geom, dipoles, gaussOrder,use_adaptive_integration)
    #ds2mm = om.DipSource2MEGMat (dipoles, squids)
    #h2mm = om.Head2MEGMat (geom, squids)
    h2em = om.Head2EEGMat (geom, electrodes)
    #gain_meg = om.GainMEG (hminv, dsm, h2mm, ds2mm)
    gain_eeg = om.GainEEG (hminv, dsm, h2em)
    gain_eeg.save(savedir)
    return gain_eeg

class LF_explore:
    def __init__(self, LF,x ,y):
        self.x = np.copy(x)
        self.y = np.copy(y)
        self.LF = LF[:]
        self.measure = np.array([])
    def addLF(self,LF):
	self.LF.append(LF)
    def sub2ind(self,i,j):
        lx = len(self.x)
        ly = len(self.y)
        return j*lx + i
    def VT(self,i,j,nv,norm):
	if nv == 8:
	    vois_i = [0, -1, -1, -1, 0, 1, 1, 1]
            vois_j = [-1, -1, 0, 1, 1, 1, 0, -1]
	else:
	    vois_i = [0, -1, 0, 1]
            vois_j = [-1, 0, 1, 0]
        lx = len(self.x)
        ly = len(self.y)
        S = 0
	if norm == 1: Sn = np.linalg.norm(self.LF[self.sub2ind(i,j)])*nv
	else: Sn = 1
        for k in range(len(vois_i)):
            ii = modulo(i + vois_i[k],lx-1)
            jj = modulo(j + vois_j[k],ly-1)	    
            S = S + np.linalg.norm(self.LF[self.sub2ind(i,j)] - self.LF[self.sub2ind(ii,jj)])/math.sqrt(math.pow(self.x[i] - self.x[ii],2) + math.pow(self.x[j] - self.x[jj],2))
        return S/Sn
    def VT2(self,i,j):
        vois_i = [0, -1, 0, 1]
        vois_j = [-1, 0, 1, 0]
        lx = len(self.x)
        ly = len(self.y)
        S = len(vois_i)*self.LF[self.sub2ind(i,j)]
        for k in range(len(vois_i)):
            ii = modulo(i + vois_i[k],lx-1)
            jj = modulo(j + vois_j[k],ly-1)
            S = S - self.LF[self.sub2ind(ii,jj)]
        return np.linalg.norm(S)
    def VTmeasure(self,mtype):
	if mtype == '4v':
	    vt = lambda i, j: self.VT(i,j,4,0)
	if mtype == '8v':
	    vt = lambda i, j: self.VT(i,j,8,0)
	if mtype == '4vn':
	    vt = lambda i, j: self.VT(i,j,4,1)
	if mtype == '8vn':
	    vt = lambda i, j: self.VT(i,j,8,1)
        for i in range(len(self.x)):
            for j in range(len(self.y)):
                self.measure = np.append(self.measure,vt(i,j))
        return self.measure
    def pointdetect(self,x,y):
	xi = np.searchsorted(self.x,x) - 1
	yi = np.searchsorted(self.y,y) - 1 
	x0 = x - self.x[xi]
	y0 = y - self.y[yi]
	cx = self.x[xi + 1] - self.x[xi]
	cy = self.y[yi + 1] - self.y[yi]
	a = cy/cx
	b = cy
	c = y0 + a*x0 - b
	if c<=0:
	    i1, j1 = xi, yi
	else:
	    i1, j1 = xi+1, yi+1
	i2, j2 = xi+1, yi
	i3, j3 = xi, yi+1
	k1, k2, k3 = self.sub2ind(i1,j1), self.sub2ind(i2,j2), self.sub2ind(i3,j3)
	A = np.array([[self.x[i1], self.x[i2], self.x[i3]],[self.y[j1], self.y[j2], self.y[j3]],[1, 1, 1]])
	b = np.array([x, y, 1])
	coef = np.linalg.solve(A,b)
	#return [[i1, i2, i3],[j1, j2, j3],coef]
	return [[k1, k2, k3], coef]
    def LFinterpol(self, x, y):
	LFnew = []
	for i in range(len(x)):
	    res = self.pointdetect(x[i],y[i])
	    idx = res[0]
	    coef = res[1]
	    LFnew.append(self.LF[idx[0]]*coef[0] + self.LF[idx[1]]*coef[1] + self.LF[idx[2]]*coef[2])
	return LFnew
