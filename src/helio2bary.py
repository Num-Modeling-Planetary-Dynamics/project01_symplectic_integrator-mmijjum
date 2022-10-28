#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:50:29 2022

@author: mmijjum

"""
import numpy as np

#This converts the velocity vector from danbys.py to barycentric

M2AU = 6.68*10**-12
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8

G = ((6.6743 * 10**-11) * M2AU**3) / (kg2Msun * sec2yr**2)  #m3/kg s2

GMsun = G #we're in units of Msun, so GMsun is just G.

GMneptune = (G) * ((1.024 * 10**26) * kg2Msun)
GMpluto = (G) * ((1.3 * 10 **22) * kg2Msun) 
GMtotal = GMneptune + GMpluto

vbvectors_neptune = []
for i in range(len(danby.vvectors_n)):
    Gmtot = GMsun + GMtotal
    vbcb = -np.sum(GMtotal * danby.vvectors_n[i], axis=0) / GMtotal + GMsun
    vbvec = danby.vvectors_n[i] + vbcb
    vbvectors_neptune.append(vbvec)

vbvectors_pluto = []
for i in range(len(danby.vvectors_p)):
    Gmtot = GMsun + GMtotal
    vbcb = -np.sum(GMtotal * danby.vvectors_p[i], axis=0) / GMtotal + GMsun
    vbvec = danby.vvectors_p[i] + vbcb
    vbvectors_pluto.append(vbvec)

