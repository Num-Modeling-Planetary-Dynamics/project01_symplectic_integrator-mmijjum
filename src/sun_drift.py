#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:44:26 2022

@author: mmijjum
"""
#This computes the linear (heliocentric) drift.

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



dt = 1/30 * danby.period

time = np.arange(0,10**5,dt)
rsunvectors = []
for i in time:
    pt = np.sum(GMtotal* helio2bary.vbvec) / GMsun
    rsunvec += pt * time[i]
    rsunvectors.append(rsunvec)
    
