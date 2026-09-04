#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 20:34:17 2026

@author: devi
"""


import math
import numpy as np
import os


filename1 = "particles.xyz"
filehandle1 = open(filename1,'r')


Ncells = 1470
Ntotal = 7000+Ncells
N = 1001
Lx = 1000.0
kd = 0.1

xp = np.zeros(Ncells,dtype=float)
zp = np.zeros(Ncells,dtype=float)
bias = np.zeros(Ncells,dtype=float)


recoc = np.zeros(N,dtype=float)

c0 = 100.0

xgrid = np.linspace(0.0,Lx,N)

skip = 10

for i in range(skip):    
    filehandle1.readline()
    filehandle1.readline()        
    for j in range(Ntotal):
        
        s = filehandle1.readline()
        split = s.split(" ")
        typ = int(split[0])
        idx = int(split[1])
        
        if typ == 1:            
            xp[idx] = float(split[2])
            zp[idx] = float(split[4])
            bias[idx] = float(split[8])
        


npar = 4000

L = np.zeros(npar,dtype=float)
time = np.zeros(npar,dtype=float)
xp_pos = xp[zp > 0.5]
L0 = (xp_pos.max() - xp_pos.min())

dt = 100/3600  # in hours

for i in range(npar):    
    filehandle1.readline()
    filehandle1.readline()        
    for j in range(Ntotal):
        
        s = filehandle1.readline()
        split = s.split(" ")
        typ = int(split[0])
        idx = int(split[1])
        
        if typ == 1:            
            xp[idx] = float(split[2])
            zp[idx] = float(split[4])
            bias[idx] = float(split[8])

    xp_pos = xp[zp > 0.5]
    L[i] = (xp_pos.max() - xp_pos.min()) / L0 if xp_pos.size else 0.0
    time[i] = dt*i
    

folder_name = os.path.basename(os.getcwd())
np.save(folder_name + ".npy", L)


