#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 12:06:20 2022

@author: mmijjum
"""
import numpy as np
import danby
import helio2bary

M2AU = 6.68*10**-12
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8

G = ((6.6743 * 10**-11) * M2AU**3) / (kg2Msun * sec2yr**2)  #m3/kg s2

GMsun = G #we're in units of Msun, so GMsun is just G.

GMneptune = (G) * ((1.024 * 10**26) * kg2Msun)
GMpluto = (G) * ((1.3 * 10 **22) * kg2Msun) 
GMtotal = GMneptune + GMpluto

mu_n = GMneptune
mu_p = GMpluto
#I skipped the HW question where we converted orbital parameters, so this is following after Dave's code
#his xv2el_one function
#but modified to theoretically work with my input parameters

rvecp = danby.rvectors_p
rvecn = danby.rvectors_n

vvecn = helio2bary.vbvectors_neptune
vvecp = helio2bary.vbvectors_pluto

rmagn = np.linalg.norm(rvecn)
rmagp = np.linalg.norm(rvecp)

vmag2n = np.vdot(vvecn,vvecn)
vmag2p = np.vdot(vvecp,vvecp)

hn = np.cross(rvecp,vvecp)
hp = np.cross(rvecn,vvecn)

hmagn = np.linalg.norm(hn)
hmagp = np.linalg.norm(hp)

rdotn = np.sign(np.vdot(rvecn,vvecn)) * np.sqrt(vmag2n - (hmagn / rmagn)**2)
rdotp = np.sign(np.vdot(rvecp,vvecp)) * np.sqrt(vmag2p - (hmagp / rmagp)**2)


an = 1.0/(2.0 / rmagn - vmag2n/mu_n)
ap = 1.0/(2.0 / rmagp - vmag2p/mu_p)

eccn = np.sqrt(1 - hmagn**2 / (mu_n * an))
eccp = np.sqrt(1 - hmagp**2 / (mu_p * ap))

incn = np.arccos(hn[2]/hmagn)
incp = np.arccos(hp[2]/hmagp)


goodincn = np.abs(incn) > np.finfo(np.float64).tiny
goodincp = np.abs(incp) > np.finfo(np.float64).tiny

sOn = np.where(goodincn,  np.sign(hn[2]) * hn[0] / (hmagn * np.sin(incn)),0.0)
sOp = np.where(goodincp,  np.sign(hp[2]) * hp[0] / (hmagp * np.sin(incp)),0.0)

cOn = np.where(goodincn, -np.sign(hn[2]) * hn[1] / (hmagn * np.sin(incn)),1.0)
cOp = np.where(goodincp, -np.sign(hp[2]) * hp[1] / (hmagp * np.sin(incp)),1.0)


Omegan = np.arctan2(sOn, cOn)
Omegap = np.arctan2(sOp, cOp)


sofn = np.where(goodincn,rvecn[2] / (rmagn * np.sin(incn)), rvecn[1]/rmagn)
cofn = np.where(goodincn,(rvecn[0] / rmagn + np.sin(Omegan) * sofn * np.cos(incn)) / np.cos(Omegan), rvecn[0]/rmagn)

sofp = np.where(goodincp,rvecp[2] / (rmagp * np.sin(incp)), rvecp[1]/rmagp)
cofp = np.where(goodincp,(rvecp[0] / rmagp + np.sin(Omegap) * sofp * np.cos(incp)) / np.cos(Omegap), rvecp[0]/rmagp)

ofn = np.arctan2(sofn,cofn)
ofp = np.arctan2(sofp,cofp)

goodeccn = eccn > np.finfo(np.float64).tiny
sfn = np.where(goodeccn,an * (1.0 - eccn**2)  * rdotn / (hmagn * eccn), 0.0)
cfn = np.where(goodeccn,(1.0 / eccn) * (an * (1.0 - eccn**2) / rmagn - 1.0), 1.0 )

goodeccp = eccp > np.finfo(np.float64).tiny
sfp = np.where(goodeccp,ap * (1.0 - eccp**2)  * rdotp / (hmagp * eccp), 0.0)
cfp = np.where(goodeccp,(1.0 / eccp) * (ap * (1.0 - eccp**2) / rmagp - 1.0), 1.0 )

fn = np.arctan2(sfn, cfn)

omegan = ofn - fn

fp = np.arctan2(sfp, cfp)

omegap = ofp - fp

omegan = np.mod(omegan, 2*np.pi)
varpin = np.mod(Omegan + omegan,2*np.pi)

omegap = np.mod(omegap, 2*np.pi)
varpip = np.mod(Omegap + omegap,2*np.pi)