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
        return i*ly + j
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
	if mtype is '4v':
	    vt = lambda i, j: self.VT(i,j,4,0)
	if mtype is '8v':
	    vt = lambda i, j: self.VT(i,j,8,0)
	if mtype is '4vn':
	    vt = lambda i, j: self.VT(i,j,4,1)
	if mtype is '8vn':
	    vt = lambda i, j: self.VT(i,j,8,1)
        for i in range(len(self.x)):
            for j in range(len(self.y)):
                self.measure = np.append(self.measure,vt(i,j))
        return self.measure
