#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:56:04 2022

@author: mmijjum
"""

#This gives you an updated velocity vector (that will get converted to barycentric in another script)
#This will also give you an updated position vector.

import numpy as np
import constants

rvec0 = constants.rvec0
vvec0 = constants.vvec0

#Get semi-major axis and eccentricity from cartesian vectors
rmag0 = np.linalg.norm(rvec0)
vmag2 = np.vdot(vvec0,vvec0)
h = np.cross(rvec0,vvec0)
hmag2 = np.dot(h,h)
a = 1/(2 / rmag0 - vmag2/mu) 
ecc = np.sqrt(1 - hmag2 / (mu*a))
P = 2 * np.pi * np.sqrt(a**3/mu) 
dt = 1/30 * P
n = np.sqrt(mu / a**3) #mean motion, kepler's third law

# Get the current value of the eccentric and mean anomalies.
# Prevent floating point exception for 0 eccentricity orbits
E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(rmag0 - a) / (a * ecc)), 0)

# If we are in the first half of the orbit (i.e. moving from periapsis to apoapsis) then the above is fine.
# However, if we are in the second half of the orbit, we have to adjust the value of E0. This is determined by
# the sign of r dot v (positive is first half, negative is second half)

if np.sign(np.vdot(rvec0, vvec0)) < 0.0:
    E0 = 2 * np.pi - E0

M0 = E0 - ecc * np.sin(E0)

# Solve Kepler's equation to advance the eccentric anomaly over a single step
M = M0 + n * dt

#Danbys:
E = []
E0 = M 
E.append(E0)
Mlist = []
while True:
    fold = E0 - ecc * np.sin(E0) - M
    fprime_old = 1 - ecc*np.cos(E0)
    fdoubleprime_old = ecc*np.sin(E0)
    ftripleprime_old = ecc*np.cos(E0)
    
    delta11 = - fold/fprime_old
    delta12 = - fold/(fprime_old + (1/2 * delta11 * fdoubleprime_old))
    delta13 = - fold / (fprime_old + (1/2 * delta12 * fdoubleprime_old) + (1/6 * delta12**2 * ftripleprime_old))

    Enew = E0 + delta13
    E0 = Enew
    M = Enew - ecc * np.sin(Enew)
    
   
    if np.absolute(Enew - E[-1]) < np.deg2rad(0.0000000000001):
        break
    E.append(Enew)
    Mlist.append(M)

E = E[-1] #Eccentric anomaly in RADIANS.
dE = E - E0

# Now implement f and g functions to advance cartesian vectors
f = a / rmag0 * (np.cos(dE) - 1.0) + 1.0
g = dt + (np.sin(dE) - dE) / n

# Advance position vector
rvec0 = np.asarray([[x],[y],[z]])
vvec0 = np.asarray([[vx], [vy], [vz]])

rvec = f * rvec0 + g * vvec0
rmag = np.linalg.norm(rvec)

fdot = -a**2 / (rmag * rmag0) * n * np.sin(dE)
gdot = a / rmag * (np.cos(dE) - 1.0) + 1.0

vvec = fdot * rvec0 + gdot * vvec0
