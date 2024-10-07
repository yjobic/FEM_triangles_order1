# -*- coding: utf-8 -*-
"""
@author: yann
"""

from FEMlib.mesh import *


def rmsError(Sol, SolEx):
    relative_error=np.sum(np.power(Sol - SolEx,2))
    return np.sqrt(relative_error/len(Sol))

def interpoleSolAtGaussPoints(ele,quad,base,Sol):
    SolInterpole=np.zeros(quad.Npts)
    for k in range(quad.Npts):
        xi= quad.coords[k][0]
        eta=quad.coords[k][1]
        #
        #[A REMPLIR]
        #
    return SolInterpole

def errorL2Elem(ele,quad,base,Sol,FuncSolEx):
    errorElem=0
    SolInterpole=interpoleSolAtGaussPoints(ele,quad,base,Sol)
    x,y=computeCoordPhyFromRef(ele, base, quad)      
    for k in range(quad.Npts):
        #here the jacobian is constant, that why we don't have to evalute it at gauss point
        #
        #[A REMPLIR]
        #
    return errorElem

def L2error(mesh, quad, Sol, FuncSolEx):
    errL2=0
    for idGroupEle in range(mesh.NbGroupEle):
        elelist=mesh.listElesType[idGroupEle]
        basis=mesh.base[idGroupEle]
        for ele in elelist:
            errL2+=errorL2Elem(ele,quad,basis,Sol,FuncSolEx)
    return np.sqrt(errL2)