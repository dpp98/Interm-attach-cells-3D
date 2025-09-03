#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 21:14:41 2025

@author: devi
"""

import math
import numpy as np
import matplotlib.pyplot as plt

filename = "particles_last.xyz"
filehandle = open(filename,'r')

filehandle2 = open("carea.xyz",'w')

Ntot = 4600
Npart = 1000

X = np.zeros(Npart,dtype=float)
Y = np.zeros(Npart,dtype=float)
Z = np.zeros(Npart,dtype=float)

zmin = -0.6

skip = 0
npar = 999

dt = 1.0

carea = np.zeros(npar,dtype=float)
time = np.zeros(npar,dtype=float)
one_area = math.pi*(0.5)**2

Lx = 59.0

total_area = Lx**2

for i in range(skip):
    
    s = filehandle.readline()
    s = filehandle.readline()
    
    for j in range(Ntot):
        
        s = filehandle.readline()

for i in range(npar):
    
    time[i] = i*dt
    
    s = filehandle.readline()
    s = filehandle.readline()
    
    for j in range(Ntot):
    
        s = filehandle.readline()
        split = s.split(" ")
    
        idx = int(split[1]) - 1
    
        if idx < Npart:

            xp = float(split[2])            
            yp = float(split[3])
            zp = float(split[4])

            if zp <= zmin + 1.0:       # Considering the cells within one cell diameter to be in-contact (Defining the wetting layer)
                carea[i] += one_area

    filehandle2.write(str(time[i])+","+str(carea[i]))
    filehandle2.write(str("\n"))
    

print(np.mean(carea),np.std(carea))

plt.plot(time,carea)
plt.show()

filehandle.close()
filehandle2.close()
