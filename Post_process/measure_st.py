# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 20:59:41 2024

@author: M K Panigrahi
"""

import math
import numpy as np
#import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

foldername = "tmpstan36"

filename = foldername+"/particles.xyz"
filehandle = open(filename,'r')

fwallhi = "fwalln36hi.dump"
fhandlehi = open(fwallhi,'r')

fwalllo = "fwalln36lo.dump"
fhandlelo = open(fwalllo,'r')

fname = foldername+"/streg.txt"
freghadl = open(fname,'a')

Npart = 4000
Ntotal = 11200
Nbc = 3600

npar = 6500
nskip = 1000


for i in range(nskip):
    
    s = filehandle.readline()
    s = filehandle.readline()
    
    for j in range(Ntotal):
        s = filehandle.readline()
    
R1 = np.zeros(npar,dtype=float)
R2 = np.zeros(2*npar,dtype=float)
R2av = np.zeros(npar,dtype=float)
R3 = np.zeros(npar,dtype=float)
sigma = np.zeros(npar,dtype=float)  

Fwalllo = np.zeros(npar,dtype=float)
Fwallhi = np.zeros(npar,dtype=float)

for i in range(npar):

    s = filehandle.readline()
    s = filehandle.readline()

    Xy0 = np.array([],dtype=float)
    Zy0 = np.array([],dtype=float)
    
    Xtop = np.array([],dtype=float)
    Ytop = np.array([],dtype=float)
    Rtop = np.array([],dtype=float)
    
    Xbot = np.array([],dtype=float)
    Ybot = np.array([],dtype=float)  
    Rbot = np.array([],dtype=float)

    for j in range(9):
        s = fhandlehi.readline()
        s = fhandlelo.readline()

    for j in range(Nbc):
        s = fhandlehi.readline()
        split = s.split(" ")
        Fwallhi[i] += float(split[7])
    
        s = fhandlelo.readline()
        split = s.split(" ")
        Fwalllo[i] -= float(split[7])

    for j in range(Ntotal):
        
        s = filehandle.readline()
        split = s.split(" ")
        idx = int(split[1])
        
        if idx < Npart:
            
            Xp = float(split[2])
            Yp = float(split[3])
            Zp = float(split[4])
            
            if Yp > -0.5 and Yp < 0.5:
                
                Xy0 = np.append(Xy0,Xp)
                Zy0 = np.append(Zy0,Zp)
                
            if Zp > 4.0 :
                Xtop = np.append(Xtop, Xp)
                Ytop = np.append(Ytop, Yp)
                Rtop = np.append(Rtop, math.sqrt(Xp**2 + Yp**2))
                
            if Zp < -4.0 :
                Xbot = np.append(Xtop, Xp)
                Ybot = np.append(Ytop, Yp)
                Rbot = np.append(Rbot, math.sqrt(Xp**2 + Yp**2))

    allPoints=np.column_stack((Xy0,Zy0))
    hullPoints = ConvexHull(allPoints)
    
    Binds = hullPoints.vertices
    
    Xb = np.array([],dtype=float)
    Zb = np.array([],dtype=float)
    
    x1 = np.array([],dtype=float)
    z1 = np.array([],dtype=float)

    x2 = np.array([],dtype=float)
    z2 = np.array([],dtype=float)
    
    for j in range(len(Binds)):
        
        Xb = np.append(Xb,Xy0[Binds[j]])
        Zb = np.append(Zb,Zy0[Binds[j]])

        if Xy0[Binds[j]] > 0.0 :
            x1 = np.append(x1,Xy0[Binds[j]])
            z1 = np.append(z1,Zy0[Binds[j]])
            
        if Xy0[Binds[j]] < 0.0 :
            x2 = np.append(x2,Xy0[Binds[j]])
            z2 = np.append(z2,Zy0[Binds[j]])

    x1mean = np.mean(x1)
    z1mean = np.mean(z1)
    
    x2mean = np.mean(x2)
    z2mean = np.mean(z2)
    
###############################    
    
    x1 = x1 - x1mean
    z1 = z1 - z1mean
    
    Suu = 0.0
    Suv = 0.0
    Svv = 0.0
    Suuu = 0.0
    Suvv = 0.0
    Svvv = 0.0
    Svuu = 0.0
    
    for j in range(len(x1)):
        
        Suu += x1[j]**2
        Svv += z1[j]**2
        Suv ++ x1[j]*z1[j]
        
        Suuu += x1[j]**3
        Svvv += z1[j]**3
        Suvv += x1[j]*z1[j]**2
        Svuu += x1[j]**2*z1[j]
        
    A = np.array([[Suu, Suv], [Suv, Svv]])
    B = np.array([0.5*(Suuu + Suvv), 0.5*(Svvv + Svuu)])
    xc = np.linalg.solve(A,B)
    
    R2[2*i] = math.sqrt(xc[0]**2 + xc[1]**2 + (Suu + Svv)/len(x1))
    
    #figure, axes = plt.subplots()
    #Drawing_uncolored_circle = plt.Circle(( xc[0]+x1mean , xc[1]+z1mean ), R2[2*i], fill = False, linestyle='--' ) 
    #axes.set_aspect( 1 )
    #axes.add_artist( Drawing_uncolored_circle )

    #plt.plot(Xy0,Zy0,'.')
    #plt.axis('equal')
    #plt.show()
    
#####################    
    
    x2 = x2 - x2mean
    z2 = z2 - z2mean
    
    Suu = 0.0
    Suv = 0.0
    Svv = 0.0
    Suuu = 0.0
    Suvv = 0.0
    Svvv = 0.0
    Svuu = 0.0
    
    for j in range(len(x2)):
        
        Suu += x2[j]**2
        Svv += z2[j]**2
        Suv ++ x2[j]*z2[j]
        
        Suuu += x2[j]**3
        Svvv += z2[j]**3
        Suvv += x2[j]*z2[j]**2
        Svuu += x2[j]**2*z2[j]
        
    A = np.array([[Suu, Suv], [Suv, Svv]])
    B = np.array([0.5*(Suuu + Suvv), 0.5*(Svvv + Svuu)])
    xc = np.linalg.solve(A,B)
        
    R2[2*i+1] = math.sqrt(xc[0]**2 + xc[1]**2 + (Suu + Svv)/len(x1))
    
    R2av[i] = 0.5*(R2[2*i] + R2[2*i+1])
    R1[i] = 0.5*(np.max(Xy0) - np.min(Xy0))
    R3[i] = 0.5*(np.max(Rtop) + np.max(Rbot))
    
    
    sigma[i] = 0.5*(Fwalllo[i] + Fwallhi[i])/(math.pi*(R3[i]**2)*(1.0/R1[i] + 1.0/R2[i]))
######################
#figure, axes = plt.subplots()
#Drawing_uncolored_circle = plt.Circle(( xc[0]+x2mean , xc[1]+z2mean ), R2av[i], fill = False, linestyle='--' )
 
#axes.set_aspect( 1 )
#axes.add_artist( Drawing_uncolored_circle )

#plt.plot(Xy0,Zy0,'.')
#plt.axis('equal')
#plt.show()

freghadl.write(foldername+"\n")
freghadl.write(str(np.mean(sigma))+" "+str(np.std(sigma))+"\n")
print(np.mean(sigma),np.std(sigma))
freghadl.close()
filehandle.close()
fhandlehi.close()
fhandlelo.close()
## mean sigma = 3.0203831773994088   std sigma = 0.5108837885582562
