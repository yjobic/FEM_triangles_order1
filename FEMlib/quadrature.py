# -*- coding: utf-8 -*-
"""
@author: yann
"""

import numpy as np


class Quadrature:
    def __init__(self,dim,ordre):
        self.dim=dim
        self.ordre=ordre
        assert dim==2, "Attention on ne gere que dim=2"
        assert ordre<5 and ordre>0,"Attention ordre doit etre en 0 et 3"
        #match instruction only valid for python > 3.10
        #Gauss quadrature, dhatt et touzot p 297
        if ordre==1:
            self.Npts=1
            a=np.double(1./3)
            self.coords=[[a,a]]
            self.wp=[np.double(0.5)]
        if ordre==2:
            self.Npts=3
            a=np.double(0.5)
            b=np.double(0)
            self.coords=[[a,a],[b,a],[a,b]]
            w=np.double(1./6)
            self.wp=[w,w,w]
        if ordre==3:
            self.Npts=4
            a=np.double(1./3)
            b=np.double(1./5)
            c=np.double(3./5)
            self.coords=[[a,a],[b,b],[c,b],[b,c]]
            w1=np.double(-27/96)
            w2=np.double(25/96)
            self.wp=[w1,w2,w2,w2]
        if ordre==4:
            self.Npts=6
            a=np.double(0.445948490915965)
            b=np.double(0.091576213509771)
            self.coords=[[a,a],[1-2*a,a],[a,1-2*a],[b,b],[1-2*b,b],[b,1-2*b]]
            w1=np.double(0.111690794839005)
            w2=np.double(0.054975871827661)
            self.wp=[w1,w1,w1,w2,w2,w2]
           
            
    def printQuad(self):
        print("dim : ",self.dim)
        print("Ordre : ",self.ordre)
        print("Nombre de points : ",self.Npts)
        print("coords : ",self.coords)
        print("poids : ",self.wp)
