import matplotlib.tri as tri
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import math


import openmeeg as om

def closest_node(node, nodes):
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)
def collinear(p0, p1, p2):
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    #print x1 * y2 - x2 * y1
    return abs(x1 * y2 - x2 * y1) > 1e-12
def supnorm(A, B):
    A = np.transpose(A)
    B = np.transpose(B)
    err = np.array([])
    for i in range(len(A)):
	err  =  np.append(err, np.linalg.norm(A[i] - B[i])/np.linalg.norm(A[i]))
    return np.amax(err)
def pinvnorm(A, B):
    A = np.linalg.pinv(A)
    B = np.linalg.pinv(B)
    return np.linalg.norm(A - B)/np.linalg.norm(A)
class TriangulationLF(tri.Triangulation):
    def __init__(self, x,y,triang = None):
        tri.Triangulation.__init__(self,x,y,triang)
        self.ver_neighbors = list([])
        self.vertice_neighbors()
        self.LF = list([])
	self.checkLF = list([])
	self.checkx = np.array([])
        self.checky = np.array([])
        self.measure = np.array([])
    def addCheckLF(self,LF,x,y):
        self.checkLF.append(LF)
        self.checkx = np.append(self.checkx,x)
        self.checky = np.append(self.checky,y)
    def addLF(self,LF):
        self.LF.append(LF)
    def add_lead_field(self, LF):
        if len(LF) != len(self.x):
            print 'Number of matrices must be same as the number of points'
            return
        self.LF = np.copy(LF)
    def combine_coef(self, x, y, ntri):
        idx = self.triangles[ntri,:]
        A = np.vstack([self.x[idx], self.y[idx], [1,1,1]])
        b = np.array([x,y,1])
        return np.linalg.solve(A,b)
    def find_triangle(self, x, y):
        t = tri.TrapezoidMapTriFinder(self)
        return t(x,y)
    def new_matrix_estim(self, x, y):
        ntri = self.find_triangle(x, y)
        idx = self.triangles[ntri,:]
	ntri = np.append(ntri,100)
	alpha = self.combine_coef(x, y, ntri[0])
	return self.LF[idx[0]]*alpha[0] + self.LF[idx[1]]*alpha[1] + self.LF[idx[2]]*alpha[2]
    def vertice_neighbors(self):
        for i in range(len(self.x)):
            neighbors = np.array([])
            for line in self.triangles:
                if i in line:
                    neighbors = np.append(neighbors,line)
            neighbors = np.unique(neighbors)
            self.ver_neighbors.append(np.delete(neighbors,np.where(neighbors==i)))
    def add_point(self, x, y,LF, methode = None):
        xx = np.append(self.x, x)
        yy = np.append(self.y, y)
        if methode is not None:
            l = len(xx)-1
            ntri = self.find_triangle(x, y)
            idx = self.triangles[ntri,:]
            triang = self.triangles
	    if collinear([xx[l],yy[l]],[xx[idx[0]],yy[idx[0]]],[xx[idx[1]],yy[idx[1]]]): triang = np.vstack((triang,[l,idx[0],idx[1]]))
            if collinear([xx[l],yy[l]],[xx[idx[2]],yy[idx[2]]],[xx[idx[1]],yy[idx[1]]]): triang = np.vstack((triang,[l,idx[1],idx[2]]))
            if collinear([xx[l],yy[l]],[xx[idx[0]],yy[idx[0]]],[xx[idx[2]],yy[idx[2]]]): triang = np.vstack((triang,[l,idx[2],idx[0]]))
            triang = np.delete(triang, ntri, axis=0)
            tri.Triangulation.__init__(self,xx,yy,triang)
        else:
            tri.Triangulation.__init__(self,xx,yy)
	self.addLF(LF)
    def check_error(self, norm):
        Err = np.array([])
        for i in range(len(self.checkx)):
            LFtmp = self.new_matrix_estim(self.checkx[i],self.checky[i])
	    if norm == 'sup': Err = np.append(Err,supnorm(self.checkLF[i],LFtmp))
	    elif norm == 'pinv': Err = np.append(Err,pinvnorm(self.checkLF[i],LFtmp))
	    else: Err = np.append(Err,np.linalg.norm(LFtmp - self.checkLF[i])/np.linalg.norm(self.checkLF[i]))
        return [np.argmax(Err), np.amax(Err), Err]
    def subtriangulation(self,x ,y):
	ntri = find_triangle(x, y)
	settri = self.neighbors[ntri,:]
    def make_grid(self, threshold,choicetype = 'refine', norm = 'L2'):
	point_set = np.transpose(np.vstack((self.checkx,self.checky)))
	i=0
	while True:
	    i+=1
	    #print i
            [midx, mval, Err] = self.check_error(norm)
            if mval < threshold:
                break
	    if choicetype == 'refine':
                self.add_point(self.checkx[midx],self.checky[midx], self.checkLF[midx])
		print mval
	    else:
		#idx = np.where(Err > threshold)
		#coeffs = Err[idx];
		#point_set_tmp = point_set[idx]
		#barycenter = np.dot(np.transpose(point_set_tmp),coeffs)/np.sum(coeffs)
		#midx = closest_node(point_set, barycenter)
		ntri = self.find_triangle(self.checkx[midx], self.checky[midx])
		ntri=np.append(ntri,100)
        	idx = self.triangles[ntri[0],:]
		coeffs = np.array([0.33,0.33,0.34])
		point_set_tmp = np.transpose(np.vstack((self.x[idx],self.y[idx])))
		barycenter = np.dot(np.transpose(point_set_tmp),coeffs)/np.sum(coeffs)
		midx = closest_node(point_set, barycenter)
		self.add_point(self.checkx[midx],self.checky[midx], self.checkLF[midx],'fixed')
	return Err
'''

    def measure_tv(self):
        for i in range(len(self.x)):
            neighbors = self.ver_neighbors[i]
            N = len(neighbors)
            Smin = 1000000
	    Smax = 0
	    for j in range(N):
                #Smin = np.minimum(Smin,np.linalg.norm(self.LF[i] - self.LF[neighbors[j]]))
	        Smax = np.maximum(Smax,np.linalg.norm(self.LF[i] - self.LF[neighbors[j]]))
	    #self.measure = np.append(self.measure, S/np.linalg.norm(self.LF[i]))
	    self.measure = np.append(self.measure, (Smax+Smax)/2)
    def smoothing(self, measure):
	new_measure = np.array([])
	for i in range(len(self.x)):
	    neighbors = self.ver_neighbors[i]
	    neighbors = [int(n) for n in neighbors]
	    allpoints = np.append(measure[neighbors],measure[i])
	    new_measure = np.append(new_measure,np.mean(allpoints))
	return new_measure
''' 
