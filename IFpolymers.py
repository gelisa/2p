#!/usr/bin/python
import numpy as np
from math import *
from scipy import special as sp

def G(x,sigma):
    return 1.0/(sigma*sqrt(2*pi))*exp(-x/(2*sigma**2))

def Delta(x,y):
    if x==y:
        return 1
    else:
        return 0

def Nkdd(l,k,df,dt):
    nkdd = 0
    for i in range(0,min(df,dt)+1):
        nkdd+=sp.binom(df,min(df,dt)-i)*sp.binom(l-df,max(0,dt-df)+i)*Delta(k,abs(dt-df)+2*i)
    
    return nkdd
    
def Pkdd(l,k,df,dt):
    return Nkdd(l,k,df,dt)/sp.binom(l,k)

def Pdd(l,df,dt,sigma):
    pdd=0
    for k in range(0,l+1):
        pdd+=Pkdd(l,k,df,dt)*G(k,sigma)
    return pdd

def CoeffB(l,b,df,sigma,rSigma,K):
    cfb=0
    for dt in range(0,l+1):
        cfb+=G(dt,rSigma)*Pdd(l,df,dt,sigma)*Pdd(l,dt,b,sigma)*K
    
    return cfb

def getA(l,sigma,rSigma,K):
    aL=[]
    for df in range(0,l+1):
        aL.append([])
        for b in range(0,l+1):
            aL[df].append(CoeffB(l,b,df,sigma,rSigma,K))
    A=np.array(aL)
    
    return A

def getPddMatrix(l,sigma):
    pM=[]
    for d1 in range(l+1):
        pM.append([])
        for d2 in range(l+1):
            pM[d1].append(Pdd(l,d1,d2,sigma))
    
    return pM

l = 10
sigma = 4.0
rSigma = 2.0
K = 1.0
b = np.array([1.0]*(l+1))
A = getA(l,sigma,rSigma,K)


x = np.linalg.solve(A,b)
print(x)
