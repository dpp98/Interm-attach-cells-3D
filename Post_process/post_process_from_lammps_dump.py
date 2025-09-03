#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:41:03 2023

@author: devi
"""

import random as rd
import numpy as np

filename1 ="ca3d_tap1_wf1_f1.dump"
file1 = open(filename1,'r')

filename2 = "ca3d_tap1_wf1_f1/particles.xyz"
file2 = open(filename2,'w')

npar = 10000
skip = 9

npart = 4600
timestep = 0.00001

x = np.zeros(npart,dtype=float)
y = np.zeros(npart,dtype=float)
z = np.zeros(npart,dtype=float)

vx = np.zeros(npart,dtype=float)
vy = np.zeros(npart,dtype=float)
vz = np.zeros(npart,dtype=float)

nid = np.zeros(npart,dtype=int)

#for i in range(2782):

#    s = file1.readline()


for i in range(npar):
    print(i)
    for j in range(skip):
        s = file1.readline()
        
    file2.write(str(npart)+"\n")
    file2.write("Timestep " + str(timestep*i)+"\n")    
        
    for j in range(npart):
        
        s = file1.readline()
        
        split = s.split(" ")
       
        typ = int(split[0])
        ind = int(split[1]) 
        x[j] = round(float(split[2]),6)
        y[j] = round(float(split[3]),6)
        z[j] = round(float(split[4]),6)
        nid[j] = int(split[6])
        vx[j] = round(float(split[8]),6)
        vy[j] = round(float(split[9]),6)
        vz[j] = round(float(split[10]),6)
        
        #file2.write(str(typ)+" "+str(ind)+" "+str(x[j])+" "+str(y[j])+" "+str(z[j])+ " " + "0.5"+"\n")
        
        file2.write(str(typ)+" "+str(ind)+" "+str(x[j])+" "+str(y[j])+" "+str(z[j])+" " + str(vx[j]) + " " + str(vy[j]) + " " + str(vz[j]) + " " + str(nid[j]) + " " + "0.5"+"\n")



file1.close()
file2.close()
