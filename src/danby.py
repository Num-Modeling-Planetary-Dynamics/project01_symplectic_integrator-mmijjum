#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:56:04 2022

@author: mmijjum
"""

import numpy as np

M2AU = 6.68*10**-12
km2AU = 6.68*10**-9
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8

G = ((6.6743 * 10**-11) * (M2AU**3)) / (kg2Msun * sec2yr**2)  #m3/kg s2

#neptune 
mu_neptune = (6835099.97 * km2AU**3) / (sec2yr**2)

#position vectors
x_n = 2.973710066866719*10**1
y_n = -3.191682991624657
z_n = -6.196085602852035*10**-1

#velocity vectors
vx_n = 3.141994243104834*10**-4
vy_n = 3.148729110427204*10**-3
vz_n = -7.207191781210733*10**-5

#pluto
mu_pluto = (869.61* km2AU**3) / (sec2yr**2)
 # from JPL horizons
x_p = 1.595092452135255*10**1
y_p = -3.070298377171116*10**1
z_p = -1.326204396461275

vx_p = 2.872759563330164*10**-3
vy_p = 7.839169969389357*10**-4
vz_p = -9.084719074113415*10**-4

#input these manuall from nssdc.gsfc.nasa.gov/planetary factsheet
#when i was calculating manually, kept getting (-) values for pluto?
a_n = 30.06
P_n= 2 * np.pi * np.sqrt(a_n**3/mu_neptune) 
dt_n = 1/30 * P_n

a_p = 39.48
P_p = 2 * np.pi * np.sqrt(a_p**3/mu_pluto) 
dt_p = 1/30 * P_p

time = np.arange

#Danbys function will calculate E, which is needed to update rvec and vvec later.
def kep(mu, dt, a, x, y, z, vx, vy, vz):
    rvec0 = [x, y, z]
    vvec0 = [vx, vy, vz]
    
    #Get semi-major axis and eccentricity from cartesian vectors
    rmag0 = np.linalg.norm(rvec0)
    #vmag2 = np.vdot(vvec0,vvec0)
    h = np.cross(rvec0,vvec0)
    hmag2 = np.dot(h,h)
    
    ecc = np.sqrt(1 - hmag2 / (mu*a))
    n = np.sqrt(mu / a**3) 
   

  # Get the current value of the eccentric and mean anomalies.
  # Prevent floating point exception for 0 eccentricity orbits
    E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(rmag0 - a) / (a * ecc)), 0)

  # If we are in the first half of the orbit (i.e. moving from periapsis to apoapsis) then the above is fine.
  # However, if we are in the second half of the orbit, we have to adjust the value of E0. This is determined by
  # the sign of r dot v (positive is first half, negative is second half)
  #took these few lines from Dave's code, modified to fit my parameters.
  
    if np.sign(np.vdot(rvec0, vvec0)) < 0.0:
        E0 = 2 * np.pi - E0

    M0 = E0 - ecc * np.sin(E0)

    E = []
    E.append(E0)
    Mlist = []
    while True:
        fold = E0 - ecc * np.sin(E0) - M0
        fprime_old = 1 - ecc*np.cos(E0)
        fdoubleprime_old = ecc*np.sin(E0)
        ftripleprime_old = ecc*np.cos(E0)
        
        delta11 = - fold/fprime_old
        delta12 = - fold/(fprime_old + (1/2 * delta11 * fdoubleprime_old))
        delta13 = - fold / (fprime_old + (1/2 * delta12 * fdoubleprime_old) + (1/6 * delta12**2 * ftripleprime_old))
    
        Enew = E0 + delta13
        E0 = Enew
        M = Enew - ecc * np.sin(Enew)
        M0 = M
        
       
        if np.absolute(Enew - E[-1]) < np.deg2rad(10**-14):
            break
        E.append(Enew)
        Mlist.append(M)

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
    #period.append(P)


    return rvec, vvec

time  = np.arange(0,10,1) #small range for now, just for testing


#this for loop will loop over different position/velocities (based on the updated position and velocities)
#from danby's method above, and theoretically give a list of new position and velocity vectors.
#doens't actually do this, coming across a divide by 0 error.

# #NEPTUNE

# rvectors_n = []
# vvectors_n = []

# for i in range(len(time)):
# rvec = danby.danbys(mu_neptune, a_n, x_n, y_n, z_n, vx_n, vy_n,vz_n)[0]
# #vvec = danby.danbys(mu_neptune, a_n, x_n, y_n, z_n, vx_n, vy_n,vz_n)[1]

# # rvectors_n.append(rvec) #this  will store the values
# # vvectors_n.append(vvec) #this will store the values

# # xnew = (np.asarray(rvec)[0])
# # ynew = (np.asarray(rvec)[1])
# # znew = (np.asarray(rvec)[2])

# # vxnew = (np.asarray(vvec)[0])
# # vynew = (np.asarray(vvec)[1]) 
# # vznew = (np.asarray(vvec)[2])

# # # x_n = xnew[0]
# # # y_n = ynew[0]
# # # z_n = znew[0]

# # # vx_n = vxnew[0]
# # # vy_n = vynew[0]
# # # vz_n = vznew[0]

# # # #PLUTO

# # # rvectors_p = []
# # # vvectors_p = []
# # # for i in range(len(time)):
# # #     rvec = danbys(G*m_pluto, x_p, y_p, z_p, vx_p, vy_p,vz_p)[0]
# # #     vvec = danbys(G*m_pluto, x_p, y_p, z_p, vx_p, vy_p,vz_p)[1]

# # #     rvectors_p.append(rvec) #this  will store the values
# # #     vvectors_p.append(vvec) #this will store the values


# # #     xnew = (np.asarray(rvec)[0])
# # #     ynew = (np.asarray(rvec)[1])
# # #     znew = (np.asarray(rvec)[2])

# # #     vxnew = (np.asarray(vvec)[0])
# # #     vynew = (np.asarray(vvec)[1]) 
# # #     vznew = (np.asarray(vvec)[2])

# # #     x_p = xnew[0]
# # #     y_p = ynew[0]
# # #     z_p = znew[0]

# # #     vx_p = vxnew[0]
# # #     vy_p = vynew[0]
# # #     vz_p = vznew[0]





                                                                              
                                                                             
    
      



   