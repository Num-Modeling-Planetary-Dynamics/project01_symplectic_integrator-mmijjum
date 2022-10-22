#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:50:29 2022

@author: mmijjum

"""
import danby
import numpy as np

#This converts the velocity vector from danbys.py to barycentric
GMsun = (6.6743 * 10**-11) * (10**30)

GMneptune = (6.6743 * 10**-11) *1.024 * 10**26
GMpluto = (6.6743 * 10**-11) * 1.3 * 10 **22
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

