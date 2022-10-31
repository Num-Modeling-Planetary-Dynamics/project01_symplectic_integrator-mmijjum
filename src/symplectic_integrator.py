#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:44:26 2022

@author: mmijjum
"""
#This computes the linear (heliocentric) drift.

import numpy as np
import matplotlib.pyplot as plt

M2AU = 6.68*10**-12
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8
km2AU = 6.68*10**-9

G = ((6.6743 * 10**-11) * M2AU**3) / (kg2Msun * sec2yr**2)  #m3/kg s2

GMsun = G #we're in units of Msun, so GMsun is just G.

GMneptune = (G) * ((1.024 * 10**26) * kg2Msun)
GMpluto = (G) * ((1.3 * 10 **22) * kg2Msun) 
GMtotal = GMneptune + GMpluto + GMsun

mu_pluto = (869.61* km2AU**3) / (sec2yr**2)
mu_neptune = (6835099.97 * km2AU**3) / (sec2yr**2)

#pluto
# from JPL horizons
xp = 1.595092452135255*10**1
yp = -3.070298377171116*10**1
zp = -1.326204396461275

vxp = 2.872759563330164*10**-3
vyp = 7.839169969389357*10**-4
vzp = -9.084719074113415*10**-4

rvec_pluto = [xp, yp, zp]
vvec_pluto = [vxp, vyp, vzp]

#neptune
xn = 2.973710066866719*10**1
yn = -3.191682991624657
zn = -6.196085602852035*10**-1

#velocity vectors
vxn = 3.141994243104834*10**-4
vyn = 3.148729110427204*10**-3
vzn = -7.207191781210733*10**-5

rvec_neptune = [xn, yn, zn]
vvec_neptune = [vxn, vyn, vzn]

vhvec= np.asarray([[0],[0],[0]])
vbcbh = -np.sum(GMtotal * vhvec) / GMtotal + GMsun
vbvecs= vhvec + vbcbh
rhvec = [0,0,0]



#Heliocentric drift: take in barycentric velocity and heliocentric position and time step, return barycentric velocity
def helio(rhvec0, vbvec, dt):
    
    
    pt = np.sum(GMtotal* vbvecs) / GMsun
    rhvec = rhvec0
    rhvec += pt * dt
    
    return rhvec

#Kick step: take in position vector of all three bodies, return sum to use for barycentric velocity kick step.
def kick(rvecn, rvecp, rvecs):
    drvecsn = rvecs - rvecn
    drvecnp = rvecn - rvecp
    irij3sn = np.linalg.norm(drvecsn, axis=1)**3
    irij3np = np.linalg.norm(drvecnp, axis=1)**3
    irij3sn = (GMsun + GMneptune) / irij3sn
    irij3np = (GMpluto + GMneptune) / irij3np

    return np.sum(irij3sn, irij3np)

vbvec += kick(rvecn, rvecp, rvecs) * dt

#Keplerian Drift: danby will update E, then using that E kep function will update position and velocity (heliocentric)
def kep(mu, rvec0, vvec0):
    #Get semi-major axis and eccentricity from cartesian vectors
    rmag0 = np.linalg.norm(rvec0)
    vmag2 = np.vdot(vvec0,vvec0)
    a = 1.0/(2.0 / rmag0 - vmag2/mu)
    P= 2 * np.pi * np.sqrt(a**3/mu) 
    dt = 1/30 * P
    n = np.sqrt(mu/a**3) 
    
    h = np.cross(rvec0,vvec0)
    hmag2 = np.dot(h,h)
    
    ecc = np.sqrt(1 - hmag2 / (mu*a))
   
    vmag2 = np.vdot(vvec0,vvec0)
    
    h = np.cross(rvec0,vvec0)
    
    hmag = np.linalg.norm(h)
    
    rdot = np.sign(np.vdot(rvec0,vvec0)) * np.sqrt(vmag2 - (hmag/ rmag0)**2)
        
    inc = np.arccos(h[2]/hmag)
    
    
    goodinc = np.abs(inc) > np.finfo(np.float64).tiny
    
    sO = np.where(goodinc,  np.sign(h[2]) * h[0] / (hmag * np.sin(inc)),0.0)
    
    cO = np.where(goodinc, -np.sign(h[2]) * h[1] / (hmag * np.sin(inc)),1.0)
    
    
    Omega = np.arctan2(sO, cO)
    
    
    sof = np.where(goodinc,rvec0[2] / (rmag0 * np.sin(inc)), rvec0[1]/rmag0)
    cof = np.where(goodinc,(rvec0[0] / rmag0 + np.sin(Omega) * sof * np.cos(inc)) / np.cos(Omega), rvec0[0]/rmag0)
    
    of = np.arctan2(sof,cof)
    
    goodecc = ecc > np.finfo(np.float64).tiny
    sf = np.where(goodecc,a * (1.0 - ecc**2)  * rdot / (hmag * ecc), 0.0)
    cf = np.where(goodecc,(1.0 / ecc) * (a * (1.0 - ecc**2) / rmag0 - 1.0), 1.0 )
    
    f = np.arctan2(sf, cf)
    
    omega = of - f
    
   
    
    omega = np.mod(omega, 2*np.pi)
    varpi = np.mod(Omega + omega,2*np.pi)
   
        
    # Get the current value of the eccentric and mean anomalies.
    E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(rmag0 - a) / (a * ecc)), 0)
    
    # If we are in the first half of the orbit (i.e. moving from periapsis to apoapsis) then the above is fine.
    # However, if we are in the second half of the orbit, we have to adjust the value of E0. This is determined by
    # the sign of r dot v (positive is first half, negative is second half)
    #took these few lines from Dave's code, modified to fit my parameters.
    
    if np.sign(np.vdot(rvec0, vvec0)) < 0.0:
        E0 = 2 * np.pi - E0
    
    M0 = E0 - ecc * np.sin(E0)

    # Solve Kepler's equation to advance the eccentric anomaly over a single step
    M = M0 + n * dt
    E = []
    dE = []
    Eold = M
    E.append(Eold)
    while True:
        fold = E0 - ecc * np.sin(E0) - M
        fprime_old = 1 - ecc*np.cos(E0)
        fdoubleprime_old = ecc*np.sin(E0)
        ftripleprime_old = ecc*np.cos(E0)
        
        delta11 = - fold/fprime_old
        delta12 = - fold/(fprime_old + (1/2 * delta11 * fdoubleprime_old))
        delta13 = - fold / (fprime_old + (1/2 * delta12 * fdoubleprime_old) + (1/6 * delta12**2 * ftripleprime_old))
        
        E = E0 + delta13
        E0 = E
        
        if np.absolute(E - E[-1]) < np.deg2rad(10**-14):
            break
        E.append(E)
        
        dE = E[-1] - E0
        dE.append(dE)
      
    f = a / rmag0 * (np.cos(dE) - 1.0) + 1.0
    g = dt + (np.sin(dE) - dE) / n
    
    # Advance position vector
    #Need to brainstorm how I can get this to be a general x/y/z instead of my planet specific x/y/z
    rvec0 = np.asarray([[x],[y],[z]])
    vvec0 = np.asarray([[vx], [vy], [vz]])
    
    rvec = f * rvec0 + g * vvec0
    rmag = np.linalg.norm(rvec)
    
    fdot = -a**2 / (rmag * rmag0) * n * np.sin(dE)
    gdot = a / rmag * (np.cos(dE) - 1.0) + 1.0
    
    vvec = fdot * rvec0 + gdot * vvec0
    
    return rvec, vvec, M, ecc, varpi

#Convert heliocentric velocity from kep step to barycentric (for helio step, next): 
def helio2bary(vhvec):
    Gmtot = GMsun + GMneptune + GMpluto
    vbcb = -np.sum(Gmtot * vhvec, axis=0) / GMtotal + GMsun
    vbvec = vhvec + vbcb
    return vbvec

#WORKFLOW: 
#First we want to run the helio kick. Inputs needed: rhvec, vbvecs (sun), and dt
#rhvec are is [0,0,0] initially (heliocentric), vbvec we can get from running helio2bary with [0,0,0]
#dt = 1/30 * period 

rhvecsun = [0,0,0]
vbvecs = helio2bary([0,0,0])
dts = 1/30 * period
dt = np.arange(0,10**5,dts)

#need to figure out a way to get this for loop to run in a way that makes sense for which planet 
#is being advanced in 'kep'. 

vbvecnlist = []
vbvecplist  = []
vbvecslist = []
rvecnlist = []
rvecplist = []
rvecslist = []
for i in dt:
    rhvec = helio(rhvecsun, vbvecs, dt)
    rhvecslist.append(rhvec)
    rvecn = kep(mu_neptune, rvec_neptune, vvec_neptune)[0]
    rvecp = kep(mu_pluto, rvec_pluto, vvec_pluto)[0]
    rvecnlist.append(rvecn)
    rvecplist.append(rvecp)
    vvecn = kep(mu_neptune, rvec_neptune, vvec_neptune)[1]
    vvecp = kep(mu_pluto, rvec_pluto, vvec_pluto)[1]
    vbvecn = helio2bary(vvecn)
    vbvecp = helio2bary(vvecp)
    vbvecs = helio2bary(vvecs)
    vbvecnlist.append(vbvecn)
    vbvecplist.append(vbvecp)
    Mn = kep(mu_neptune, rvec_neptune, vvec_neptune)[2]
    Mp = kep(mu_pluto, rvec_pluto, vvec_pluto)[2]
    En = kep(mu_neptune, rvec_neptune, vvec_neptune)[3]
    Ep = kep(mu_pluto, rvec_pluto, vvec_pluto)[3]
    
    #now update?

#ok at some point the stuff above is going to give me the updated position/velocities
#assuming that's true, calculate energy

#Energy
ken = []
pen = []

kep = []
pep = []

kes = []
pes = []
#Neptune
for i in range(len(vbvecnlist)):
    ke = 1/2*(GMneptune) * (vbvecnlist[i]**2)
    pe = -GMneptune*rvecn[i]
    ken.append(ke)
    pen.append(pe)

#Pluto
for i in range(len(vbvecplist)):
    ke = 1/2*(GMpluto) * (vbvecplist[i]**2)
    pe = -GMpluto*rvecp[i]
    kep.append(ke)
    pep.append(pe)
    
#Sun
for i in range(len(vbvecslist)):
    ke = 1/2*(GMsun) * (vbvecslist[i]**2) #need barycentric solar velocity
    pe = -GMsun*rhvec[i]
    kes.append(ke)
    pes.append(pe)


total_e = ken + pen + kep + pep + kes + pes
    
    
varpin = kep(mu_neptune, rvec_neptune, vvec_neptune)[4]
varpip = kep(mu_pluto, rvec_pluto, vvec_pluto)[4]


lamn = Mn + varpin
lamp = Mp + varpip

omega = 3 * lamp - 2*lamn - varpin

dt = 1/30 * period

time = np.arange(0,10**5,dt)

plt.plot(time, dE/E0)
plt.plot(time, omega)


