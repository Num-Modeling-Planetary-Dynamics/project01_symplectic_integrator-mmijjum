#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 10:55:33 2022

@author: mmijjum
"""
import danby
import numpy as np
# Now implement f and g functions to advance cartesian vectors
M2AU = 6.68*10**-12
km2AU = 6.68*10**-9
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8

G = ((6.6743 * 10**-11) * (M2AU**3)) / (kg2Msun * sec2yr**2)  #m3/kg s2

#neptune 
mu_neptune = (6835099.97 * km2AU**3) / (sec2yr**2)

#position vectors
x_n = 2.973710066866719*10**1
y_n = -3.191682991624657
z_n = -6.196085602852035*10**-1

#velocity vectors
vx_n = 3.141994243104834*10**-4
vy_n = 3.148729110427204*10**-3
vz_n = -7.207191781210733*10**-5

#pluto
mu_pluto = (869.61* km2AU**3) / (sec2yr**2)
 # from JPL horizons
x_p = 1.595092452135255*10**1
y_p = -3.070298377171116*10**1
z_p = -1.326204396461275

vx_p = 2.872759563330164*10**-3
vy_p = 7.839169969389357*10**-4
vz_p = -9.084719074113415*10**-4

#input these manuall from nssdc.gsfc.nasa.gov/planetary factsheet
#when i was calculating manually, kept getting (-) values for pluto?
a_n = 30.06
P_n= 2 * np.pi * np.sqrt(a_n**3/mu_neptune) 
dt_n = 1/30 * P_n

a_p = 39.48
P_p = 2 * np.pi * np.sqrt(a_p**3/mu_pluto) 
dt_p = 1/30 * P_p

