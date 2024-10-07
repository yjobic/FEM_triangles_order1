# -*- coding: utf-8 -*-
"""
@author: yann
"""

import sys

#local import
from FEMlib.Triplets import *
from FEMlib.GeomFactors import *

class Basis:
    def __init__():
        self.dim=0
        self.Nphi=0  #On a Nphi fonction de forme pour ndof
    
class LagrangeT1(Basis):
    name='Base T1 lagrange'
    
    def __init__(self,dim,ndof):
        self.dim=dim
        self.Nphi=ndof #On a Nphi fonction de forme pour ndof
        assert dim==2, "Attention on ne gere que dim=2"
        
    def Phi1(xi,eta):
        return 1-eta-xi
    def Phi2(xi,eta):
        return xi
    def Phi3(xi,eta):
        return eta
    FctForm=[Phi1,Phi2,Phi3]
    
    #On passe aux derivées des fonctions de formes dans l'espace de reference
    def gradPhi(self, i) :
        if(i == 0):
            return np.array([-1,-1])
        elif (i == 1):
            return np.array([1,0])
        return np.array([0,1])

    #On definit ici la construction des matrices elementaire
    #a part car comme les elements sont lineaire, il y a de grosses
    #simplifications lors d'integration numerique
    def constructMassElem(self, ele, quad):
        matrix = np.zeros((self.Nphi,self.Nphi))
        for k in range(quad.Npts):
            xi=quad.coords[k][0]
            eta=quad.coords[k][1]
            for i in range(self.Nphi):
                for j in range(self.Nphi):
                    matrix[i,j] += quad.wp[k]*ele.detJ \
                                  *self.FctForm[i](xi,eta) \
                                  *self.FctForm[j](xi,eta)
        return matrix

    def constructRigidElem(self, ele, quad):
        matrix = np.zeros((self.Nphi,self.Nphi))
        Bp = computeInvJ(computeJac(ele).transpose())
        BpTBp=np.matmul(Bp.transpose(),Bp)
        coeff = ele.vol
        for i in range(self.Nphi):
            gphi_i=np.matmul(BpTBp,self.gradPhi(i))
            for j in range(self.Nphi):
                val=coeff*np.matmul(self.gradPhi(j).transpose(),gphi_i)
                matrix[i,j]=val
        return matrix
    
    
    def constructVectSolElem(self, mesh, ele, base, quad, fonc):
        vect = np.zeros(len(ele.p))                
        x,y=computeCoordPhyFromRef(ele, base, quad)      
        for k in range(quad.Npts):
            xi =quad.coords[k][0]
            eta=quad.coords[k][1]
            for i in range(self.Nphi):
                phi1=self.FctForm[i](xi,eta) 
                vect[i]+=quad.wp[k] * ele.detJ * phi1 * fonc(x[k],y[k])
        return vect














