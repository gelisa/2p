#!/usr/bin/python
import numpy as np
from math import *
from scipy import special as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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

def dFdt(l,m,sigma,ki,kf,I,F):
    dFdt=[]
    for d in range(l+1):
        dFddt=m-F[d]
        for i in range(l+1):
            dFddt+=ki*I[i]*Pdd(l,i,d,sigma)
        dFdt.append(dFddt)
    
    return dFdt

def dIdt(l,m,sigma,rSigma,ki,kf,I,F):
    dIdt=[]
    for d in range(l+1):
        dIddt=m-I[d]
        for i in range(l+1):
            dIddt+=kf*F[i]*G(i,rSigma)*Pdd(l,i,d,sigma)
        dIdt.append(dIddt)
    
    return dIdt

def dydt(l,m,sigma,rSigma,ki,kf,y):
    I=y[0:l+1]
    F=y[l+1:]
    return dIdt(l,m,sigma,rSigma,ki,kf,I,F)+dFdt(l,m,sigma,ki,kf,I,F)


l=15
m=4
sigma=3.0
rSigma=0.5
ki=2.0
kf=2.0
t=np.linspace(0, 10., 100)
def f(y,t):
    global m
    global l
    global sigma
    global rSigma
    global ki
    global kf
    return dydt(l,m,sigma,rSigma,ki,kf,y)
    
y0=[10.0]*(2*(l+1))


soln = odeint(f,y0, t)
I=[]
F=[]
for i in range(l+1):
    I.append(soln[:,i])
for i in range(l+1,2*(l+1)):
    F.append(soln[:,i])
    

plt.figure()
for i in range(l+1):
    plt.plot(t,I[i],label='I['+str(i)+']')
plt.legend(loc=0)
plt.show()