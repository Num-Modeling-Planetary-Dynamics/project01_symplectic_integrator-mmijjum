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

rsunvec = [0,0,0]

GMsun = (6.6743 * 10**-11) * (10**30)
GMneptune = (6.6743 * 10**-11) *1.024 * 10**26
GMpluto = (6.6743 * 10**-11) * 1.3 * 10 **22
GMtotal = GMneptune + GMpluto


dt = 1/30 * danby.period

time = np.arange(0,10**5,dt)
rsunvectors = []
for i in time:
    pt = np.sum(GMtotal* helio2bary.vbvec) / GMsun
    rsunvec += pt * time[i]
    rsunvectors.append(rsunvec)
    
