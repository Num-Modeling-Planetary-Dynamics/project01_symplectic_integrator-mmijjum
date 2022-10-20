#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 17:04:15 2022

@author: mmijjum
"""
import helio2bary
import numpy as np
import constants
import danby

#This is supposed to calculate the potential and kinetic energy of the system
#it clearly doesn't do that, though

#but this is needed to calculate the Hamiltonian of the system (total energy, which is what I'm supposed to plot)
vbcb = helio2bary.vbec
vbmag2 = (vbcb[0]**2+vbcb[1]**2+vbcb[2]**2)

irh = 1.0 / np.linalg.norm(constants.rvec, axis=1)
ke = 0.5 * (constants.GMcb * np.vdot(vbcb, vbcb) + np.sum(constants.Gmass * vbmag2)) / constants.GC


drvec = danby.rvec[i + 1:, :] - danby.rvec[i, :]
irij = np.linalg.norm(drvec, axis=1)
irij = constants.Gm[i + 1:] * constants.Gm[i] / irij

n = constants.rvec.shape[0]
pe = (-constants.GMcb * np.sum(constants.Gmass * irh) - np.sum([pe for i in range(n - 1)])) / constants.GC

total_E = ke + pe


