#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:14:51 2024

@author: devi
"""

import math
import numpy as np
#import matplotlib.pyplot as plt

filename1 = "particles.xyz"
filehandle1 = open(filename1,'r')


Ncells = 1470
Ntotal = 7000+1470
N = 1001
npar = 1500
kd = 0.1

dt = 100/60  # in minutes

xp = np.zeros(Ncells,dtype=float)
zp = np.zeros(Ncells,dtype=float)
bias = np.zeros(Ncells,dtype=float)

time = np.zeros(npar,dtype=float)
av_bias = np.zeros(npar,dtype=float)
av_xpos = np.zeros(npar,dtype=float)

recoc = np.zeros(N,dtype=float)

c0 = 100.0

skip = 500

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

for i in range(npar):
    
    filehandle1.readline()
    filehandle1.readline()    
    
    xmax = 0.0
    xmin = 1000.0
    
    drop_x = []
    drop_bias = []
    
    for j in range(Ntotal):
        
        s = filehandle1.readline()
        split = s.split(" ")
        typ = int(split[0])
        idx = int(split[1])
        
        if typ == 1:
            
            xp[idx] = float(split[2])
            zp[idx] = float(split[4])
            bias[idx] = float(split[8])
            
            if zp[idx] > 0.5:    # considering a thershold of 15um for the drop
                
                drop_x.append(xp[idx])
                drop_bias.append(bias[idx])
            	

    av_bias[i] = np.mean(drop_bias)
     
    av_xpos[i] = np.mean(drop_x)  

drop_vel =  np.zeros(npar-100,dtype=float)

for i in range(npar-100):
    
    drop_vel[i] = 10.0* (av_xpos[i+100] - av_xpos[i])/(100*dt) # units of microns per min
 
print(np.mean(av_bias),np.mean(drop_vel),np.std(drop_vel))
    
filehandle1.close()
