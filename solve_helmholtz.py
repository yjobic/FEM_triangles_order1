# -*- coding: utf-8 -*-
"""
@author: yann Jobic
"""

from scipy.sparse.linalg import spsolve
import matplotlib.pyplot  as plt
import numpy as np
import time

#local import
from FEMlib.mesh import *
from FEMlib.assembleProblem import  *
from FEMlib.quadrature import  *
from FEMlib.basis import *
from FEMlib.errors import *


def g(x,y):
    #
    #[A REMPLIR]
    #
def f(x,y):
    #
    #[A REMPLIR]
    #

def printTemps(temps):
    if temps>3600: 
        return '{:2.3f}'.format(temps/3600)+' h'
    elif temps>60: 
        return '{:2.3f}'.format(temps/60)+' min'
    else: 
        return '{:2.3f}'.format(temps)+' s'

plotMesh=1
outPlot=0
if  outPlot :
    get_ipython().run_line_magic('matplotlib', 'qt')
else :
    get_ipython().run_line_magic('matplotlib', 'inline')

start = time.time()

#load the mesh
#ORDER 1 mandatory
MeshFileName="Meshes/square_tri_o1.msh"
mesh = Mesh()
startloc = time.time()
mesh.GmshToMesh(MeshFileName)
endloc = time.time(); elapsed = endloc - startloc
print('Temps GmshToMesh : '+ printTemps(elapsed))

#On positionne les espaces d'approximations internes
mesh.updateBase()

#On choisit la quadrature utilisee pour l'integration
quad = Quadrature(2,3)

#t is the maxtrix stored in the COO format (for COOrdiantes)
startloc = time.time()
t=assembleMatrix(mesh, 0, quad, 1, 0)
endloc = time.time(); elapsed = endloc - startloc
print('Temps assembleMatrix : ' + printTemps(elapsed))

startloc = time.time()
B=assembleVecForceElem(mesh, 0, quad, f)
endloc = time.time(); elapsed = endloc - startloc
print('Temps assembleVecForceElem : ' + printTemps(elapsed))

#Dirichlet boundary condition applied
Dirichlet(mesh,t,g,B)

#We solve the linear system
A=coo_matrix(t.data).tocsr()
X = spsolve(A, B) # solve avec scipy

#Exact solution at nodes
SolEx=np.zeros(mesh.Npts)
for i in range(mesh.Npts):
    SolEx[i]=g(mesh.points[i].coord[0],mesh.points[i].coord[1])


# Compute errors
rmsE = rmsError(X, SolEx)
print(f'RMS error : {rmsE}')
errL2 = L2error(mesh,quad,X,g)
print(f'L2 error : {errL2}')
print('Nombre total d elements : ',mesh.getTotEle())


end = time.time()
elapsed = end - start
print('Temps d\'exécution : ' + printTemps(elapsed))

x= [pt.coord[0] for pt in mesh.points]
y= [pt.coord[1] for pt in mesh.points]
connectivity=[]
for tri in mesh.listElesType[0]:
    connectivity.append([ p.id for p in tri.p])
plt.tricontourf(x, y, connectivity, X, 12)
if plotMesh:
    plt.triplot(x, y, connectivity, 'k-', lw=1.0)
plt.colorbar()
plt.show()

