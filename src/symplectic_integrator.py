#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:44:26 2022

@author: mmijjum
"""
#This computes the linear (heliocentric) drift.

import numpy as np
import os
import matplotlib.pyplot as plt

os.chdir(os.path.dirname(os.path.abspath(__file__)))
directory = os.path.dirname(__file__)

M2AU = 6.68*10**-12
kg2Msun = 5.05*10**-31
sec2yr = 3.17*10**-8
km2AU = 6.68*10**-9

GC =  6.67428 * 10**(-11) #m3/kg s2
AU2M = 1.49597871e+11
GMSunSI = 1.3271244e+20
MSun = 1.98840987e+30
JD2S = 86400

MSun = 1.0,
MNeptune = 19412.24
MPluto = 1.35e8
Gmass = []
msun_over_mpl = [1, 19412.24, 1.35e8] #sun, neptune, pluto
GMSunAUD = GMSunSI * JD2S**2 / AU2M**3

for i in range(len(msun_over_mpl)):
    Gmasstemp = GMSunAUD/msun_over_mpl[i]
    Gmass.append(Gmasstemp)

Gmtot = GMSunAUD + np.sum(Gmass)


GMNeptuneAUD = GC*MNeptune * JD2S**2 / AU2M**3
GMPlutoAUD = GC*MPluto * JD2S**2 / AU2M**3
GMtotal = GMNeptuneAUD + GMPlutoAUD + GMSunAUD

#FORMAT: [SUN, NEPTUNE, PLUTO]

mu = [1, GMSunAUD * (1.0 + 1.0/MNeptune), GMSunAUD * (1.0 + 1.0/MPluto)]

# from JPL horizons
x = [0, 2.973710066866719*10**1, 1.595092452135255*10**1]
y = [0, -3.191682991624657, -3.070298377171116*10**1]
z = [0, -6.196085602852035*10**-1, -1.326204396461275]

vx = [0, 3.141994243104834*10**-4, 2.872759563330164*10**-3]
vy = [0, 3.148729110427204*10**-3, 7.839169969389357*10**-4]
vz = [0,  -7.207191781210733*10**-5, -9.084719074113415*10**-4]

#vbcbh = -np.sum(GMtotal * vhvec) / GMtotal + GMSunAUD
#vbvecs = vhvec + vbcbh
rhvec = [0,0,0]

    
def xv2el(mu, x,y,z,vx,vy,vz):
    rvec = ([x,y,z])
    vvec = ([vx,vy,vz])
    rmag = np.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    vmag2 = np.vdot(vvec, vvec)
    
    h = np.cross(rvec, vvec)
    hmag = np.linalg.norm(h)
    hx = h[0]
    hy = h[1]
    hz = h[2]
    
    rdot = np.sign(np.vdot(rvec,vvec)) * np.sqrt(vmag2 - (hmag / rmag)**2)
    
    a = (1.0/(2.0 / rmag - vmag2/mu))
    
    P = 2 * np.pi * np.sqrt(a**3/mu)
    
    e = (np.sqrt(1 - (hmag**2/(mu*a))))
    
    inc = np.rad2deg(np.arccos(hz/hmag))
    
    goodinc = np.abs(np.deg2rad(inc)) > np.finfo(np.float64).tiny
    sO = np.where(goodinc,  np.sign(hz) * hx / (hmag * np.sin((np.deg2rad(inc)))),0.0)
    cO = np.where(goodinc, -np.sign(hz) * hy / (hmag * np.sin((np.deg2rad(inc)))),1.0)
    Omega = (np.rad2deg(np.arctan2(sO, cO)))
    goodecc = e > np.finfo(np.float64).tiny
    sf = np.where(goodecc,a * (1.0 - e**2)  * rdot / (hmag * e), 0.0)
    cf = np.where(goodecc,(1.0 / e) * (a * (1.0 - e**2) / rmag - 1.0), 1.0 )
    sof = np.where(goodinc,rvec[2] / (rmag * np.sin(np.deg2rad(inc))), rvec[1]/rmag)
    cof = np.where(goodinc,(rvec[0] / rmag + np.sin(np.deg2rad(Omega)) * sof * np.cos(np.deg2rad(inc))) / np.cos(np.deg2rad(Omega)), rvec[0]/rmag)
    of = np.arctan2(sof,cof)
    
    f = (np.rad2deg(np.arctan2(sf, cf)))
    
    omega = (np.rad2deg(np.mod(of - np.deg2rad(f), 2*np.pi)))

    varpi = (np.rad2deg(np.mod(np.deg2rad(Omega) + np.deg2rad(omega),2*np.pi)))

    E = np.where(e > np.finfo(np.float64).tiny, np.arccos(-(rmag - a) / (a * e)), 0)
   
    if np.sign(np.vdot(rvec, vvec)) < 0.0:
        
        E = 2 * np.pi - E
    
    M = (np.rad2deg(E - e * np.sin(E)))


    lam = (np.rad2deg(np.mod(np.deg2rad(M) + np.deg2rad(varpi),2*np.pi)))


    return rvec, vvec, a, e, inc, f, omega, E, M, lam, P, varpi


#Heliocentric Drift
#Convert heliocentric velocity from kep step to barycentric (for helio step, next): 
def helio2bary(vhvec):
    vbcb = -np.sum(Gmtot * np.array(vhvec), axis=0) / Gmtot + GMSunAUD
    vbvec = vhvec + vbcb
    return vbvec,vbcb #vbcb = central body, vbvec = other body


def helio(Gmass,vbvec,rhvec):
    pt = np.sum(Gmass * vbvec) / GMSunAUD
    rhvec += pt * dt
    return rhvec


def kick(Gmass,rhvec_list):
   
    drvecsun_neptune = rhvec_list[0] - rhvec_list[1]
    drvecneptune_pluto = rhvec_list[1] - rhvec_list[2]
    drvecsun_pluto = rhvec_list[0] - rhvec_list[2]
    
    drvec = [drvecsun_neptune,drvecneptune_pluto,drvecsun_pluto]
    
    irij3 = np.linalg.norm(drvec, axis = 1)**3
    xs = []
    for i in range(len(drvec)):
        irij3[i] = 1
        irij3 = Gmass[i] / irij3
        x =  np.sum(np.array(drvec).T * irij3, axis = 1)
        xs.append(x)
    
    return xs

def accel(Gmass,rhvec_list,vbvec):
        
    vbvec += np.array(kick(Gmass,rhvec_list)) * dt
    
    return vbvec

def calc_energy(mu, vbvec_list, rhvec):
    #Energy
    for i in range(len(vbvec_list)):
        ke = 1/2 * mu * vbvec_list[i]**2
        pe = -mu * rhvec[i]
        total_e = ke + pe
    return total_e
    

if __name__ == '__main__':
   
    time = np.arange(0,10**5)
    E = []
    o_plot = []
    for j in range(len(time)):
        elements = []
        for i in range(len(mu)):
            temp = xv2el(mu[i],x[i], y[i], z[i], vx[i], vy[i], vz[i])
            elements.append(temp)
            
        sun_el = elements[0]
        neptune_el = elements[1]
        pluto_el = elements[2]
        
        
        rvec_list = []
        rvec_list.append(sun_el[0])
        rvec_list.append(neptune_el[0])
        rvec_list.append(pluto_el[0])
        
        
        vvec_list = []
        vvec_list.append(sun_el[1])
        vvec_list.append(neptune_el[1])
        vvec_list.append(pluto_el[1])
        
        Plist = []
        Plist.append(elements[1][10])
        Plist.append(elements[2][10])
        
        Pfinal = np.min(Plist)
        
        dt = Pfinal/30
        
        vbvec = helio2bary(vvec_list)[0]
        vbcb = helio2bary(vvec_list)[1]
        rhvec_list = []
        rhvec_list.append(helio(Gmass,vbvec,rvec_list)[0]) #sun
        rhvec_list.append(helio(Gmass,vbvec,rvec_list)[1]) #neptune
        rhvec_list.append(helio(Gmass,vbvec,rvec_list)[2]) #pluto
        vbvec_list = []
        vbvec_list.append(accel(Gmass,rhvec_list,vbvec))
        for i in range(len(mu)):
            E.append(calc_energy(mu[i],vbvec_list,rhvec_list))
       
        o = 3 * elements[2][9] - 2 * elements[1][9] - elements[2][11]
        o_plot.append(o)
       
        x = [rhvec_list[0][0], rhvec_list[1][0], rhvec_list[2][0]]
        y = [rhvec_list[0][1], rhvec_list[1][1], rhvec_list[2][1]]
        z = [rhvec_list[0][2], rhvec_list[1][2], rhvec_list[2][2]]
        
    plt.plot(E)
    plt.plot(o, time)
        
        


    
