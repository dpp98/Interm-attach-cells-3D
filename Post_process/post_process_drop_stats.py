#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 15:20:14 2026

@author: devi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 20:34:17 2026

@author: devi
"""


import math
import numpy as np
import os
import matplotlib.pyplot as plt


filename1 = "particles.xyz"
filehandle1 = open(filename1,'r')

filename2 = "chem_conc.xyz"
filehandle2 = open(filename2,'r')

Ncells = 1470
Ntotal = 7000+Ncells
N = 1001
Lx = 1000.0
kd = 2.0
dx = Lx / (N-1)

xp = np.zeros(Ncells,dtype=float)
zp = np.zeros(Ncells,dtype=float)
bias = np.zeros(Ncells,dtype=float)
occpc = np.zeros(Ncells,dtype=float)

c = np.zeros(N,dtype=float)
dc = np.zeros(N,dtype=float)
recoc = np.zeros(N,dtype=float)

c0 = 100.0

xgrid = np.linspace(0.0,Lx,N)

skip = 1

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
            
    filehandle2.readline()
    filehandle2.readline()
    
    s = filehandle2.readline()
    split = s.split(' ')
    for j in range(N):
        c[j] = float(split[j])
        
    s = filehandle2.readline()
    split = s.split(' ')
    for j in range(N):
        dc[j] = float(split[j])
        


npar = 4000

Ldrop = np.zeros(npar,dtype=float)
L = np.zeros(npar,dtype=float)
avgrocc = np.zeros(npar,dtype=float)
stdrocc = np.zeros(npar,dtype=float)

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
            
    filehandle2.readline()
    filehandle2.readline()
    
    s = filehandle2.readline()
    split = s.split(' ')
    for j in range(N):
        c[j] = float(split[j])
        
    s = filehandle2.readline()
    split = s.split(' ')
    for j in range(N):
        dc[j] = float(split[j])
        
    for j in range(Ncells):
        
        ind = round(xp[j]/dx)
        occpc[j] = c[ind] / (c[ind] + kd) 

    avgrocc[i] = np.mean(occpc)
    stdrocc[i] = np.std(occpc)


    xp_pos = xp[zp > 0.5]
    
    Ldrop[i] = 10.0*(xp_pos.max() - xp_pos.min()) if xp_pos.size else 0.0
    L[i] = (xp_pos.max() - xp_pos.min()) / L0 if xp_pos.size else 0.0
    time[i] = dt*i
    

# -----------------------------
# Plot aggregate size and receptor occupancy
# -----------------------------

fig, ax1 = plt.subplots(figsize=(6.0, 4.0))

# Aggregate size -- left y-axis
ax1.plot(
    time,
    Ldrop,
    color='midnightblue',
    linestyle=':',
    linewidth=1.0
)

ax1.set_xlabel('Time (Hours)', fontsize=14)
ax1.set_ylabel('Aggregate size ($\u03bcm$)', fontsize=14, color='midnightblue')

ax1.tick_params(
    axis='y',
    labelcolor='midnightblue',
    labelsize=11
)
ax1.tick_params(axis='x', labelsize=11)

# Match the approximate limits in the example figure
ax1.set_xlim(0, 110)
ax1.set_ylim(170, 205)


# Average receptor occupancy -- right y-axis
ax2 = ax1.twinx()

ax2.plot(
    time,
    avgrocc,
    color='green',
    linewidth=1.8
)

ax2.set_ylabel(
    'Average receptor occupancy',
    fontsize=14,
    color='green'
)

ax2.tick_params(
    axis='y',
    labelcolor='green',
    labelsize=11
)

ax2.set_ylim(0.8, 1.0)


# General formatting
ax1.tick_params(direction='out')
ax2.tick_params(direction='out')

plt.tight_layout()

plt.savefig('drop_stats.svg')

plt.show()


