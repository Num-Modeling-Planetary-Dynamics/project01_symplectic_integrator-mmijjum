#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:46:52 2022

@author: mmijjum
"""
#These are some constants and initialized vectors we need

x = 2.973710066866719*10**1
y = -3.191682991624657
z = -6.196085602852035*10**-1
rx = x
ry = y 
rz = z

vx = 3.141994243104834*10**-4
vy = 3.148729110427204*10**-3
vz = -7.207191781210733*10**-5

rvec0 = [rx, ry, rz]
vvec0 = [vx, vy, vz]

GC = 6.6743e-11 #m3/kg s2
Mcb = 1.988e30 
GMcb = GC*Mcb
MSun_over_Mpl = 19412.24 #Neptune
#MSun_over_Mpl = 1.35e8 #Pluto

Gmass = GMcb/ MSun_over_Mpl
Mbody = 1e26
mu = GC*(Mcb + Mbody) 
