#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 23:12:19 2025

@author: devi
"""

import math
import numpy as np
#import matplotlib.pyplot as plt

filename1 = "particles.xyz"
filehandle1 = open(filename1,'r')

filename2 = "xvout.txt"
filehandle2 = open(filename2,'w')

filename3 = "chem_conc.xyz"
filehandle3 = open(filename3,'r')



Ncells = 1470
Ntotal = 7000+1470
npar = 3930      # number of timesnaps considered
N = 1001
Lx = 1000.0
kd = 2.0
dx = Lx / (N-1)

xp = np.zeros(Ncells,dtype=float)
yp = np.zeros(Ncells,dtype=float)
zp = np.zeros(Ncells,dtype=float)
bias = np.zeros(Ncells,dtype=float)
avgrocc = np.zeros(npar,dtype=float)
stdrocc = np.zeros(npar,dtype=float)
xc = np.zeros(npar,dtype=float)
yc = np.zeros(npar,dtype=float)
zc = np.zeros(npar,dtype=float)
avb = np.zeros(npar,dtype=float)
occpc = np.zeros(Ncells,dtype=float)
bias = np.zeros(Ncells,dtype=float)
c = np.zeros(N,dtype=float)
dc = np.zeros(N,dtype=float)
Ldrop = np.zeros(npar,dtype=float)
time = np.zeros(npar,dtype=float)

dt = 100/60.0   # in minutes based on the frequency of output


for i in range(npar):
    
    filehandle1.readline()
    filehandle1.readline()    
    
    xmax = 0.0
    xmin = 1000.0    
    
    for j in range(Ntotal):
        
        s = filehandle1.readline()
        split = s.split(" ")
        typ = int(split[0])
        idx = int(split[1])
        
        if typ == 1:
            
            xp[idx-1] = float(split[2])
            zp[idx-1] = float(split[4])
            bias[idx-1] = float(split[8])
            
            if zp[idx-1] > 0.5:
                if xp[idx-1] > xmax:
                    xmax = xp[idx-1]
                if xp[idx-1] < xmin:
                    xmin = xp[idx-1]
            
    
    nsam = 0
    
    for j in range(Ncells):
        
        if zp[j] > 0.0:     # ignore the wetting layer
            xc[i] += xp[j]
            yc[i] += yp[j]
            zc[i] += zp[j]
            avb[i] += bias[j]
            nsam += 1
            
    xc[i] = xc[i] / nsam
    yc[i] = yc[i] / nsam
    zc[i] = zc[i] / nsam
    avb[i] = avb[i] / nsam
    
    Ldrop[i] = xmax - xmin
    time[i] = i*dt
    
    filehandle3.readline()
    filehandle3.readline()
    
    s = filehandle3.readline()
    split = s.split(' ')
    for j in range(N):
        c[j] = float(split[j])
        
    s = filehandle3.readline()
    split = s.split(' ')
    for j in range(N):
        dc[j] = float(split[j])
    
    for j in range(Ncells):
        
        ind = round(xp[j]/dx)
        occpc[j] = c[ind] / (c[ind] + kd) 

    avgrocc[i] = np.mean(occpc)
    stdrocc[i] = np.std(occpc)
    

for i in range(npar):
    
    filehandle2.write(str(xc[i])+","+str(yc[i])+","+str(zc[i])+","+str(avb[i])+","+str(Ldrop[i])+","+str(avgrocc[i])+","+str(stdrocc[i]))
    filehandle2.write("\n")


bias_steady = avb[npar-1]

for i in range(npar):
    
    if abs(avb[i]-bias_steady) < 0.02:
        start = i
        break
        
avx = []
step = 50
for i in range(start+step,npar,step):
    
    avx.append(10.0*(xc[i]-xc[i-step])/(step*dt))    # Here we adjust the units of velocity to be microns/min
    

aver_vel = np.mean(avx)
std_vel = np.std(avx)

print(start,npar)
print(len(avx))
print(bias_steady,aver_vel,std_vel,avgrocc[npar-1],stdrocc[npar-1])

#plt.plot(time,xc,'r')
#plt.show()

filehandle1.close()
filehandle2.close()
filehandle3.close()
