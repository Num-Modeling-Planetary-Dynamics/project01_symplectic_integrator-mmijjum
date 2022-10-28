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


#input these manually from nssdc.gsfc.nasa.gov/planetary factsheet
#when i was calculating manually, kept getting (-) values for pluto?
a_n = 30.06
P_n= 2 * np.pi * np.sqrt(a_n**3/mu_neptune) 
dt_n = 1/30 * P_n


time = np.arange(0, 10**5, dt_n)

def kep(x, y, z, vx, vy, vz):
    rvec0 = [x, y, z]
    vvec0 = [vx, vy, vz]
    
    #Get semi-major axis and eccentricity from cartesian vectors
    rmag0 = np.linalg.norm(rvec0)
    #vmag2 = np.vdot(vvec0,vvec0)
    n = np.sqrt(mu_neptune / a_n**3) 
    
    
    h = np.cross(rvec0,vvec0)
    hmag2 = np.dot(h,h)
    
    ecc = np.sqrt(1 - hmag2 / (mu_neptune*a_n))
    
    
    # Get the current value of the eccentric and mean anomalies.
    E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(rmag0 - a_n) / (a_n * ecc)), 0)
    
    # If we are in the first half of the orbit (i.e. moving from periapsis to apoapsis) then the above is fine.
    # However, if we are in the second half of the orbit, we have to adjust the value of E0. This is determined by
    # the sign of r dot v (positive is first half, negative is second half)
    #took these few lines from Dave's code, modified to fit my parameters.
    
    if np.sign(np.vdot(rvec0, vvec0)) < 0.0:
        E0 = 2 * np.pi - E0
    
    M0 = E0 - ecc * np.sin(E0)
    
    E = []
    Mlist = []
    E.append(E0)
    
    
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
    
    f = a_n / rmag0 * (np.cos(E) - 1.0) + 1.0
    g = dt_n + (np.sin(E) - E) / n
    
    # Advance position vector
    rvec0 = np.asarray([[x_n],[y_n],[z_n]])
    vvec0 = np.asarray([[vx_n], [vy_n], [vz_n]])
    
    rvec = f * rvec0 + g * vvec0
    rmag = np.linalg.norm(rvec)
    
    fdot = -a_n**2 / (rmag * rmag0) * n * np.sin(E)
    gdot = a_n / rmag * (np.cos(E) - 1.0) + 1.0
    
    vvec = fdot * rvec0 + gdot * vvec0
    
    return rvec, vvec

#NEPTUNE

rvectors_n = []
vvectors_n = []

##This loop is supposed to go into the above function, update the x/y/z/vx/vy/vz, and re-run
##Then the rvectors_n and vvectors_n is supposed to save it so I can use it in other scripts
##However, it is getting mad, encountering the following error: 
#/Users/mmijjum/Desktop/third_year/Midterm/project01_symplectic_integrator-mmijjum/src/kep_neptune.py:54: 
# RuntimeWarning: invalid value encountered in sqrt ecc = np.sqrt(1 - hmag2 / (mu_neptune*a_n))
    
for i in time:
    rvec = kep(x_n, y_n, z_n, vx_n, vy_n,vz_n)[0]
    vvec = kep(x_n, y_n, z_n, vx_n, vy_n,vz_n)[1]
    
    rvectors_n.append(rvec) #this  will store the values
    vvectors_n.append(vvec) #this will store the values
    
    xnew = (np.asarray(rvec)[0])
    ynew = (np.asarray(rvec)[1])
    znew = (np.asarray(rvec)[2])
    
    vxnew = (np.asarray(vvec)[0])
    vynew = (np.asarray(vvec)[1]) 
    vznew = (np.asarray(vvec)[2])
    
    x_n = xnew[0]
    y_n = ynew[0]
    z_n = znew[0]

    vx_n = vxnew[0]
    vy_n = vynew[0]
    vz_n = vznew[0]

