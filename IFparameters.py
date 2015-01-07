#!/home/eliza/anaconda3/bin/python3.4
import os
L = [15, 30]
M = [1.0, 3.0, 10.0]
SIGMA = [0.3, 1.0, 3.0, 10.0]
RSIGMA = [0.3, 1.0, 3.0, 10.0]
KI = [1.0, 3.0, 10.0]
KF = [1.0, 3.0, 10.0]
File=''
def writesh(l,m,sigma,rSigma,ki,kf):
    filename= 's'+str(l)+'-'+str(m)+'-'+str(sigma)+'-'+str(rSigma)+'-'+str(ki)+'-'+str(kf)+'.sh'
    global File
    File = open(filename, 'w')
    def getout(str):
        global File
        File.write(str + '\n')
    
    getout('#!/bin/bash')
    getout('#$ -S /bin/bash')
    getout('#$ -N twoP')
    getout('#$ -cwd')
    getout('#$ -q ' + 'cpu_short')
    getout('#$ -pe openmpi 1')
    getout('#$ -P kenprj')
    getout('/cavern/eliza/origins/twoPolymers/IFode.py '+str(l)+' '+str(m)+' '+str(sigma)+' '+str(rSigma)+' '+str(ki)+' '+str(kf))

    File.close()
    print('.sh written')
    return filename

for l in L:
    for m in M:
        for sigma in SIGMA:
            for rSigma in RSIGMA:
                for ki in KI:
                    for kf in KF:
                        scriptName=writesh(l,m,sigma,rSigma,ki,kf)
                        os.system('chmod ug+x '+scriptName)
                        os.system('./'+scriptName)
                        