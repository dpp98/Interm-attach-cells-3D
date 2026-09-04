#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 09:23:51 2025

@author: devi
"""

import math
import numpy as np
#import matplotlib.pyplot as plt

filename = "Dt.xyz"
filehandle = open(filename,'r')

nlines = 8000
oosix = 1.0/6.0

time = np.zeros(nlines,dtype=float)
msd = np.zeros(nlines,dtype=float)

for i in range(nlines):
    
    s = filehandle.readline()
    
    split = s.split(" ")
    
    time[i] = float(split[0])
    msd[i] = float(split[1])


nskip = 3000
dt = 100.0 # in s
deff = np.zeros(nlines-nskip,dtype=float)
stride = 500

for i in range(nlines-nskip):
    
    deff[i] = 100 * oosix * (msd[nskip+i]-msd[nskip+i-stride]) / (stride*dt) # micron^2/s


print(np.mean(deff),np.std(deff))


filehandle.close()
