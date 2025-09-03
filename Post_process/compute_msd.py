#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 20:55:51 2025

@author: devi
"""

import math
import numpy as np
#import matplotlib.pyplot as plt

filename = "particles.xyz"
filehandle = open(filename,'r')

fileout = open("Dt.xyz",'w')

Ncells = 1327
npar = 7700

time = np.zeros(npar-1,dtype=float)
xp = np.zeros((npar,Ncells),dtype=float)
yp = np.zeros((npar,Ncells),dtype=float)
zp = np.zeros((npar,Ncells),dtype=float)

xpuw = np.zeros((npar,Ncells),dtype=float)
ypuw = np.zeros((npar,Ncells),dtype=float)
zpuw = np.zeros((npar,Ncells),dtype=float)

dt = 0.001
Lx = 11.0

for i in range(npar):
    print("reading ",i)
    s = filehandle.readline()  
    s = filehandle.readline()
    for j in range(Ncells):
        
        s = filehandle.readline()
        split = s.split(" ")
        idx = int(split[1])
        x = float(split[2])
        y = float(split[3])
        z = float(split[4])
        
        xp[i][idx-1] = x
        yp[i][idx-1] = y
        zp[i][idx-1] = z
            
xpuw[0] = xp[0]
ypuw[0] = yp[0]
zpuw[0] = zp[0]

for i in range(1,npar):
    print("correcting ",i) 
    for j in range(Ncells):
        
        deltax = xp[i][j] - xp[i-1][j]
        deltay = yp[i][j] - yp[i-1][j]
        deltaz = zp[i][j] - zp[i-1][j]
                
        if deltax > 0.5*Lx:
            deltax -= Lx
        elif deltax < -0.5*Lx:
            deltax += Lx

        if deltay > 0.5*Lx:
            deltay -= Lx
        elif deltay < -0.5*Lx:
            deltay += Lx
            
        if deltaz > 0.5*Lx:
            deltaz -= Lx
        elif deltaz < -0.5*Lx:
            deltaz += Lx
            
        xpuw[i][j] = xpuw[i-1][j] + deltax
        ypuw[i][j] = ypuw[i-1][j] + deltay
        zpuw[i][j] = zpuw[i-1][j] + deltaz
    
msd = np.zeros(npar-1,dtype=float)
    
for i in range(1,npar):
    for j in range(Ncells):
        
        time[i-1] = i*dt
        msd[i-1] += math.sqrt((xpuw[i][j] - xpuw[0][j])**2 + (ypuw[i][j] - ypuw[0][j])**2 + (zpuw[i][j] - zpuw[0][j])**2)
        
    msd[i-1] /= Ncells
    
    fileout.write(str(i)+" "+str(msd[i-1]))
    fileout.write("\n")
    


#plt.plot(time,Dt)
#plt.xlim([3.0,5.5])
#plt.ylim([50.0,70.0])
#plt.show()    
    
filehandle.close()
fileout.close()
    
    
    
