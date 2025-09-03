#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 16:02:28 2025

@author: devi
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull  # <-- Added for convex hull

filename = "particles_last.xyz"
filehandle = open(filename,'r')

my_dpi=100
Ntot = 4600
Npart = 1000

X = np.zeros(Npart,dtype=float)
Y = np.zeros(Npart,dtype=float)
Z = np.zeros(Npart,dtype=float)


zmin = -0.6 

npar = 1000 # number of snapshots

dt = 1.0

area = np.zeros(npar,dtype=float)
perm = np.zeros(npar,dtype=float)
shape = np.zeros(npar,dtype=float)
Req = np.zeros(npar,dtype=float)


dmax = 6.0

for i in range(npar):
    
    
    s = filehandle.readline()
    s = filehandle.readline()
    
    xcluster = []
    ycluster = []
    
    xcluster2 = []
    ycluster2 = []
    
    # get the cluster based on z cordinate
    for j in range(Ntot):
    
        s = filehandle.readline()
        split = s.split(" ")
    
        idx = int(split[1]) - 1
    
        if idx < Npart:

            xp = float(split[2])            
            yp = float(split[3])
            zp = float(split[4])
            
            if zp > zmin + 1.0:
                
                xcluster = np.append(xcluster,xp)
                ycluster = np.append(ycluster,yp)
        
        
    # clean the cluster
    xcenter = np.mean(np.asarray(xcluster))
    ycenter = np.mean(np.asarray(ycluster))
    for j in range(len(xcluster)):
            
        dist = math.sqrt((xcluster[j]-xcenter)**2 + (ycluster[j]-ycenter)**2)
            
        if dist < dmax:
            xcluster2 = np.append(xcluster2,xcluster[j])
            ycluster2 = np.append(ycluster2,ycluster[j])

    
    ## Get the centroid of the cluster ##
    xcenter = np.mean(np.asarray(xcluster2))
    ycenter = np.mean(np.asarray(ycluster2))

    # ---- Convex Hull Section ----
    if len(xcluster2) >= 3:  # Need at least 3 points for a hull
        points = np.column_stack((xcluster2, ycluster2))
        hull = ConvexHull(points)
        
                    
    # Compute area using shoelace formula
    hull_points = points[hull.vertices]
    x = hull_points[:, 0]
    y = hull_points[:, 1]
    area[i] = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    # Perimeter (length of convex hull edges)
    perm[i] = np.sum(
        np.sqrt(
            np.sum(np.diff(np.vstack([hull_points, hull_points[0]]), axis=0)**2, axis=1)
            )
        )

    
    Req[i] = math.sqrt(area[i]/math.pi)
    
    
    # Points of the convex hull
    hull_points = points[hull.vertices]
    xbdary = hull_points[:, 0]
    ybdary = hull_points[:, 1]
    dist = np.asarray(xbdary)    
    
    xcenter = np.mean(np.asarray(xbdary))
    ycenter = np.mean(np.asarray(ybdary))
    
    
    ## Get the centroid of the cluster ##
    xcenter_clus = np.mean(np.asarray(xcluster2))
    ycenter_clus = np.mean(np.asarray(ycluster2))
    
    for j in range(len(xbdary)):
        
        dist[j] = math.sqrt((xbdary[j]-xcenter)**2 + (ybdary[j]-ycenter)**2)
    
    
    shape[i] = np.std(dist)
    
    if i % 10 == 0:
        fig = plt.figure(figsize=(6, 6), dpi=my_dpi)  # 600 x 600 pixels
        plt.scatter(xcluster2, ycluster2, c='green')
        plt.xlim([xcenter_clus-6.0,xcenter_clus+6.0])
        plt.ylim([ycenter_clus-6.0,ycenter_clus+6.0])


        for simplex in hull.simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], 'r--')

        plt.axis('off')  # Turn off axes
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Remove all margins

        if i < 100:
            filename = "shape" + "0" + str(int(i/10)+1) + ".png"
        else:
            filename = "shape" + str(int(i/10)+1) + ".png"

        plt.savefig(filename)  # Save exactly as configured (600x600 px, no cropping)
        plt.show()
        plt.close()





print(np.mean(shape)/np.mean(Req), np.std(shape)/np.mean(Req))
