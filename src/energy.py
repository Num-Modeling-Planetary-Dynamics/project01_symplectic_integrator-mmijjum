#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:07:35 2022

@author: mmijjum
"""

import numpy as np
import danby
import helio2bary
import sun_drift

vbvn = helio2bary.vbvectors_neptune
vbvp = helio2bary.vbvectors_pluto

rn = danby.rvectors_n
rp = danby.rvectors_p
rsun = sun_drift.rsunvectors

M2AU = 6.68*10**-12
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8

G = ((6.6743 * 10**-11) * M2AU**3) / (kg2Msun * sec2yr**2)  #m3/kg s2

GMsun = G #we're in units of Msun, so GMsun is just G.

GMneptune = (G) * ((1.024 * 10**26) * kg2Msun)
GMpluto = (G) * ((1.3 * 10 **22) * kg2Msun) 
GMtotal = GMneptune + GMpluto


#Energy
ken = []
pen = []

kep = []
pep = []

kes = []
pes = []
#Neptune
for i in range(len(rn)):
    ke = 1/2*(GMneptune) * (vbvn[i]**2)
    pe = -GMneptune*rn[i]
    ken.append(ke)
    pen.append(pe)

#Pluto
for i in range(len(rp)):
    ke = 1/2*(GMpluto) * (vbvp[i]**2)
    pe = -GMpluto*rp[i]
    kep.append(ke)
    pep.append(pe)
    
#Sun
for i in range(len(rsun)):
    ke = 1/2*(GMsun) * (helio2bary.vbvs[i]**2) #need barycentric solar velocity
    pe = -GMsun*rsun[i]
    kes.append(ke)
    pes.append(pe)


total_e = ken + pen + kep + pep + kes + pes