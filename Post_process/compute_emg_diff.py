#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 09:23:51 2025

@author: devi
"""

import math
import numpy as np
import matplotlib.pyplot as plt

filename = "Dt.xyz"
filehandle = open(filename,'r')

nlines = 9999
oosix = 1.0/6.0

time = np.zeros(nlines,dtype=float)
msd = np.zeros(nlines,dtype=float)

for i in range(nlines):
    
    s = filehandle.readline()
    
    split = s.split(" ")
    
    time[i] = float(split[0])*(100/3600)
    msd[i] = float(split[1])*(100/1000000)


nskip = 5000       # Skip initial few frames to ensure that the measurement is taken at steady state
dt = 100.0/3600.0
deff = np.zeros(nlines-nskip,dtype=float)
stride = 500

for i in range(nlines-nskip):
    
    deff[i] = (1000000/60) * oosix * (msd[nskip+i]-msd[nskip+i-stride]) / (stride*dt) # micron^2/min


print(np.mean(deff),np.std(deff))


filehandle.close()
