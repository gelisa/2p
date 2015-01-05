#!/usr/bin/python
import numpy as np
from math import *
from scipy import special as sp


def G(x,sigma):
    '''Gaussian distribution'''
    return 1.0/(sigma*sqrt(2*pi))*exp(-x/(2*sigma**2))

def Delta(x,y):
    '''cronecker delta'''
    if x==y:
        return 1
    else:
        return 0

def Nkdd(l,k,df,dt):
    '''number of sequences which are on the dt shell in distance k away from df shell'''
    nkdd = 0
    for i in range(0,min(df,dt)+1):
        nkdd+=sp.binom(df,min(df,dt)-i)*sp.binom(l-df,max(0,dt-df)+i)*Delta(k,abs(dt-df)+2*i)
    
    return nkdd
    
def Pkdd(l,k,df,dt):
    '''probability to make a jump of length k from df shell into dt shell
    '''
    return Nkdd(l,k,df,dt)/sp.binom(l,k)

def Pdd(l,df,dt,sigma):
    '''probability to make a jump from df to dt'''
    pdd=0
    for k in range(0,l+1):
        pdd+=Pkdd(l,k,df,dt)*G(k,sigma)
    return pdd

def CoeffB(l,b,df,sigma,rSigma,ki,kf):
    ''' an element of the matrix of coefficients for an ode, produced from the probabilities above '''
    cfb=0
    for dt in range(0,l+1):
        cfb+=G(dt,rSigma)*Pdd(l,df,dt,sigma)*Pdd(l,dt,b,sigma)*ki*kf
    
    return cfb

def getA(l,sigma,rSigma,ki,kf):
    ''' returns matrix mentioned above '''
    aL=[]
    for df in range(0,l+1):
        aL.append([])
        for b in range(0,l+1):
            if not b==df:
                aL[df].append(CoeffB(l,b,df,sigma,rSigma,ki,kf))
            else:
                aL[df].append(CoeffB(l,b,df,sigma,rSigma,ki,kf)-1)
    A=np.array(aL)
    
    return A

def getPddMatrix(l,sigma):
    ''' returns the matrix of Pdd's mentioned above '''
    pM=[]
    for d1 in range(l+1):
        pM.append([])
        for d2 in range(l+1):
            pM[d1].append(Pdd(l,d1,d2,sigma))
    
    return pM

def getB(l,m,kf,rSigma):
    B=[]
    for d in range(l+1):
        Bd=-m
        for i in range(l+1):
            Bd+=m*kf*G(i,rSigma)*Pdd(l,i,d,sigma)
        B.append(Bd)
        
    
    return np.array(B) 

l=5
sigma=1.0
rSigma=0.5
ki=2.0
kf=3.0
m=2.0

A=getA(l,sigma,rSigma,ki,kf)
B=getB(l,m,kf,rSigma)

    
x = np.linalg.solve(A, B)    
print(x)



