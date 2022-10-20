#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:50:29 2022

@author: mmijjum
"""
#This converts the velocity vector from danbys.py to barycentric

import danby
import constants
import numpy as np

vhvec = danby.vvec
#vh2vb: 
vbcb = -np.sum(constants.Gmass * vhvec, axis=0) / constants.Gmtot
vbvec = vhvec + vbcb

# vb2vh:
vbcb = -np.sum(constants.Gmass * vbvec, axis=0) / constants.GMcb
vhvec = vbvec - vbcb
