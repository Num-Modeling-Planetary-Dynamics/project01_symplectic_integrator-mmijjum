#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:44:26 2022

@author: mmijjum
"""
#This computes the linear (heliocentric) drift.

import helio2bary
import constants
import numpy as np

vbvec = helio2bary.vbvec
Gmass = constants.Gmass
GMcb = constants.GMcb
dt = constants.dt

rhvec = constants.rvec0
#Linear drift

pt = np.sum(Gmass * vbvec) / GMcb
rhvec += pt * dt

#this is the updated position vector for the Sun.