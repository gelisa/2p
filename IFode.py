#!/usr/bin/python
import numpy as np
from math import *
from scipy import special as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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