# -*- coding: utf-8 -*-
"""
@author: yann
"""

import numpy as np
from FEMlib.mesh import *


def computeJac(ele):
    Jac=np.array([[np.double(0),np.double(0)],[np.double(0),np.double(0)]])
    #i -> ligne : (dx,dy)
    for i in range(2):
        #j -> colonne : (dxi,deta)
        for j in range(2):
            #ici ele triangulaire, forme analytique
            Jac[i][j]=ele.p[j+1].coord[i]-ele.p[0].coord[i]
    return Jac

def computeCoordPhyFromRef(ele, base, quad):
    x=np.zeros(quad.Npts)
    y=np.zeros(quad.Npts)
    #On cherche la position des points de gauss dans l'espace physique
    for k in range(quad.Npts):
        xi= quad.coords[k][0]
        eta=quad.coords[k][1]
        for i in range(len(ele.p)):
            x[k]+=base.FctForm[i](xi,eta)*ele.p[i].coord[0]
            y[k]+=base.FctForm[i](xi,eta)*ele.p[i].coord[1]
    return x,y
            
def computeDetJ(self,ele):
    ele.detJ = np.linalg.det(ele.jac)

def computeInvJ(Jac):
    return np.linalg.inv(Jac)
        