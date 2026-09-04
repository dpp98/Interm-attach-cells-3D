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

fileout = open("speed.xyz",'w')

Ncells = 1327
npar = 8000

time = np.zeros(npar-1,dtype=float)
xp = np.zeros((npar,Ncells),dtype=float)
yp = np.zeros((npar,Ncells),dtype=float)
zp = np.zeros((npar,Ncells),dtype=float)

xpuw = np.zeros((npar,Ncells),dtype=float)
ypuw = np.zeros((npar,Ncells),dtype=float)
zpuw = np.zeros((npar,Ncells),dtype=float)

dt = 100.0  # in s 
Lx = 11.0

for i in range(npar):
    #print("reading ",i)
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
    #print("correcting ",i) 
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
    
Dt = np.zeros(npar-1,dtype=float)


skip = 3000
gap = 500

speed = np.zeros(npar-skip,dtype=float)    

for i in range(skip,npar):
    for j in range(Ncells):
        
        
        delx = xpuw[i][j] - xpuw[i-gap][j]
        dely = ypuw[i][j] - ypuw[i-gap][j]
        delz = zpuw[i][j] - zpuw[i-gap][j]
        
        time[i-1] = i*dt
        
        speed[i-skip] += 60*10.0*math.sqrt(delx**2 + dely**2 + delz**2) / (gap*dt) # in um/min
        
    speed[i-skip] /= Ncells
    
    fileout.write(str(i)+" "+str(speed[i-skip]))
    fileout.write("\n")
    

print(np.mean(speed),np.std(speed))
       
filehandle.close()
fileout.close()
    
    
    
