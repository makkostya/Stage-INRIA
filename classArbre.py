import numpy as np
from classTriangulationLF import *
def deviser_tout(s):
    if not s.sommets:
        s.deviser()
    else:
        for som in s.sommets:
            deviser_tout(som)
def unique_points(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
def error_calcul(x,y,LF,a,seuil,norm = 'L2'):
    xx, yy = np.meshgrid(x, y)
    xx = xx.flatten()
    yy = yy.flatten()
    go = True
    j = 1
    z = 1
    while go:
	Err_rectangle = np.array([])
        go = False
	j+=1
        for i in range(len(xx)):
            p = np.array([xx[i], yy[i]])
            inf = a.acces_sommet(p,'interpol')
            ps = inf[0]
	    coeff = inf[1]
	    LFest = 0
            for k in range(3):
                LFest += get_LF(x,y,LF,ps[k])*coeff[k]
	    if norm == 'sup': err = supnorm(get_LF(x,y,LF,p),LFest)
	    elif norm == 'pinv': err = pinvnorm(get_LF(x,y,LF,p),LFest)
	    else: err = np.linalg.norm(LFest - get_LF(x,y,LF,p))/np.linalg.norm(get_LF(x,y,LF,p))
	    Err_rectangle = np.append(Err_rectangle, err)
            if err >= seuil:
                go = True
                a.acces_sommet(p,'deviser')
    return Err_rectangle
def get_LF(x,y, LF, p):
    ilast = len(x)-1
    xi = np.searchsorted(x,p[0],side='right') - 1
    d1 = p[0] - x[xi]
    d2 = x[min(xi+1,ilast)] - p[0]
    if d1>=d2: xi = min(xi+1,ilast)
    yi = np.searchsorted(y,p[1],side='right') - 1
    d1 = p[1] - y[yi]
    d2 = y[min(yi+1,ilast)] - p[1]
    if d1>=d2: yi = min(yi+1,ilast)
    i = xi + yi*len(x)
    return LF[i]
class GridArbre:
    def __init__(self, v1, v2, code = 'r'):
        self.v1, self.v2 = v1, v2
        self.sommets = []
        self.code = code
    def deviser(self):
	if len(self.sommets)>0:
	    print 'Alarme!!'
	    return
        iy = [0, 0, 1, 1]
        ix = [0, 1, 0, 1]
        w = self.v2[0] - self.v1[0]
        h = self.v2[1] - self.v1[1]
        for i in range(4):
            v1tmp = self.v1 + np.array([ix[i]*w/2.0, iy[i]*h/2.0])
            v2tmp = v1tmp + np.array([w/2.0, h/2.0])
            sommet = GridArbre(v1tmp,v2tmp,str(i))
            self.sommets.append(sommet)
    def acces_sommet(self,x,tache = None):
	if len(self.sommets)>0:
            for sommet in self.sommets:
                if (sommet.v1[0]-1e-15)<=x[0]<=(sommet.v2[0]+1e-15) and (sommet.v1[1]-1e-15)<=x[1]<=(sommet.v2[1]+1e-15):
                    if tache == 'adresse': return self.code + sommet.acces_sommet(x,'adresse')
                    if tache == 'interpol': return sommet.acces_sommet(x,'interpol')
                    if tache == 'deviser': return sommet.acces_sommet(x,'deviser')
	else:
            if tache == 'adresse':
		print self.code
		return self.code
            if tache == 'interpol': 
		tmp = self.interpolation(x)
		return tmp
            if tache == 'deviser':
		self.deviser()
		return
    def get_rectangles(self):
        tmp = []
        if len(self.sommets):
            for sommet in self.sommets:
                tmp+=sommet.get_rectangles()
            return tmp
        w = self.v2[0] - self.v1[0]
        h = self.v2[1] - self.v1[1]
        return [[self.v1, w, h]]
    def interpolation(self,x):
        w = self.v2[0] - self.v1[0]
        h = self.v2[1] - self.v1[1]
        x1 = self.v1 + np.array([0, h])
        x2 = self.v1 + np.array([w, 0])
        a = (x1[1] - x2[1])/(x1[0]-x2[0])
        b = x1[1] - a*x1[0]
        if x[1] - a*x[0] <= b : t = self.v1
        else: t = self.v2
        A = np.transpose(np.vstack([t, x1, x2]))
        A = np.vstack([A,[1,1,1]])
        B = np.append(x,1)
        return [np.vstack([t, x1, x2]), np.linalg.solve(A,B)]
    def get_points(self):
	tmp = np.empty((0,2),float)
        if len(self.sommets):
            for sommet in self.sommets:
                tmp = np.vstack((tmp, sommet.get_points()))
            return tmp
	#a = np.array([ [self.v1[0], self.v1[1]], [self.v2[0], self.v1[1]] , [self.v1[0], self.v2[1]] , [self.v2[0], self.v2[1]] ])
	#print np.vstack((a, a))
        return np.array([ [self.v1[0], self.v1[1]], [self.v2[0], self.v1[1]] , [self.v1[0], self.v2[1]] , [self.v2[0], self.v2[1]] ])
