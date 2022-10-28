#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 10:55:33 2022

@author: mmijjum
"""
import danby
import numpy as np
# Now implement f and g functions to advance cartesian vectors

def kep(a, E, rmag0, n, dt):
    f = a / rmag0 * (np.cos(E) - 1.0) + 1.0
    g = dt + (np.sin(E) - E) / n
    
    # Advance position vector
    rvec0 = np.asarray([[x],[y],[z]])
    vvec0 = np.asarray([[vx], [vy], [vz]])
    
    rvec = f * rvec0 + g * vvec0
    rmag = np.linalg.norm(rvec)
    
    fdot = -a**2 / (rmag * rmag0) * n * np.sin(E)
    gdot = a / rmag * (np.cos(E) - 1.0) + 1.0
    
    vvec = fdot * rvec0 + gdot * vvec0
    period.append(P)

    
    return rvec, vvec

time  = np.arange(0,10,1) #small range for now, just for testing


#this for loop will loop over different position/velocities (based on the updated position and velocities)
#from danby's method above, and theoretically give a list of new position and velocity vectors.
#doens't actually do this, coming across a divide by 0 error.

#NEPTUNE

rvectors_n = []
vvectors_n = []

for i in range(len(time)):
    rvec = danbys(G*m_neptune, x_n, y_n, z_n, vx_n, vy_n,vz_n)[0]
    vvec = danbys(G*m_neptune, x_n, y_n, z_n, vx_n, vy_n,vz_n)[1]
    
    rvectors_n.append(rvec) #this  will store the values
    vvectors_n.append(vvec) #this will store the values
    
    xnew = (np.asarray(rvec)[0])
    ynew = (np.asarray(rvec)[1])
    znew = (np.asarray(rvec)[2])
    
    vxnew = (np.asarray(vvec)[0])
    vynew = (np.asarray(vvec)[1]) 
    vznew = (np.asarray(vvec)[2])

    # x_n = xnew[0]
    # y_n = ynew[0]
    # z_n = znew[0]
    
    # vx_n = vxnew[0]
    # vy_n = vynew[0]
    # vz_n = vznew[0]
    
# # #PLUTO

# # rvectors_p = []
# # vvectors_p = []
# # for i in range(len(time)):
# #     rvec = danbys(G*m_pluto, x_p, y_p, z_p, vx_p, vy_p,vz_p)[0]
# #     vvec = danbys(G*m_pluto, x_p, y_p, z_p, vx_p, vy_p,vz_p)[1]
    
# #     rvectors_p.append(rvec) #this  will store the values
# #     vvectors_p.append(vvec) #this will store the values

    
# #     xnew = (np.asarray(rvec)[0])
# #     ynew = (np.asarray(rvec)[1])
# #     znew = (np.asarray(rvec)[2])
    
# #     vxnew = (np.asarray(vvec)[0])
# #     vynew = (np.asarray(vvec)[1]) 
# #     vznew = (np.asarray(vvec)[2])

# #     x_p = xnew[0]
# #     y_p = ynew[0]
# #     z_p = znew[0]
    
# #     vx_p = vxnew[0]
# #     vy_p = vynew[0]
# #     vz_p = vznew[0]
    
    


    
                                                                              
                                                                             
    
    