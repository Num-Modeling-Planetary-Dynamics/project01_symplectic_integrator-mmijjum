#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 11:56:40 2022

@author: mmijjum
"""

import numpy as np
import energy
import danby
import orbital_elements
import matplotlib.pyplot as plt

#Plot dE/E0:
    
E = energy.total_e 

dE = []
E0 = []

lam_neptune = danby.M +orbital_elements.varpin

lam_pluto = danby.M +orbital_elements.varpip

omega = 3 * lam_pluto - 2*lam_neptune - orbital_elements.varpip


for i in range(len(E)):
    Eold = E[i-1]
    E = E[i]
    
    deltaE = Eold - E
    dE.append(deltaE)
    E0.append(Eold)

dt = 1/30 * danby.period

time = np.arange(0,10**5,dt)

plt.plot(time, dE/E0)
plt.plot(time, omega)