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

GMsun = (6.6743 * 10**-11) * (10**30)
GMneptune = (6.6743 * 10**-11) *1.024 * 10**26
GMpluto = (6.6743 * 10**-11) * 1.3 * 10 **22
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