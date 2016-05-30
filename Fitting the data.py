# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:12:48 2016

Fitting model to data

@author: findlaywebb
"""

import numpy as np
import Model as md
import fittingFunctions as fit
import pylab as plt
import pandas as pd
from scipy.optimize import minimize
from os import getcwd

## Load the data

subRoute = getcwd()
transitData = fit.load(subRoute + r"/Qatar 2 transit 10:04:16/Object data/reReduced/",r"Measurements.xls")

time = transitData["J.D.-2400000"]
timeR = time - time[0]
relFlux = transitData["rel_flux_T1"]
relFluxErr = transitData["rel_flux_err_T1"]

def model(t, p_rad, p_mass, p_t0, scale, offset):
    """
    t       -   is the time values
    p_rad   -   is the radius of the planet
    p_mass  -   is the mass of the planet
    p_t0    -   sets when the first transit happens
    scale   -   scales the lightcuvres relative flux
    offset  -   Uniformally increases or decreases the relative flux
    """    
    
    s_mass = 0.74 #Mass of parent star in stellar masses
    s_rad = 0.713 #Radius of parent star in stellar radii
    p_period = 1.33712 #Period of the planet as it can't be calculated from one transit    
    
    lc = md.lightcurve([s_mass,s_rad],[p_rad,p_mass,p_period,p_t0],t)
    lc=(lc-1)*scale + 1 - offset
    
    return lc

##Initial guesses [radius, mass, t0, scale, offset]
x0_t0 = 0.0497191844603
x0 = [0.117512,0.002375,x0_t0,10,-9.7]

m0 = minimize(fit.ModelChiSquared, x0, method = 'Nelder-Mead', args=(model,[timeR,relFlux,relFluxErr],), options={'maxiter' : 1e10,'disp': True})

x1 = m0['x']

for i in np.arange(0,3,1):
    m = minimize(fit.ModelChiSquared, x1, method = 'Nelder-Mead', args=(model,[timeR,relFlux,relFluxErr],), options={'maxiter' : 1e10,'disp': True})

def printResults(boo=True):
    if boo:
        print 'Success = ',m['success']
        print 'Chi-squared = ',m['fun']

        DoF = len(relFluxErr)-2
        print 'Reduced chi-squared = ',m['fun']/DoF

        print 'Planet radius = ',m['x'][0]
        print 'Planet mass = ',m['x'][1]
        print 'Planet t0 = ',m['x'][2]
        print 'Scale = ',m['x'][3]
        print 'Offset = ',m['x'][4]
        
printResults(True)

n=25000
Params = fit.Bootstrap([timeR,relFlux,relFluxErr],m['x'],model,n)
print 'Radius:'
plt.hist(Params[0],bins=15, color = 'darkorange')
plt.ylabel("Count")
plt.xlabel("Radius /Solar radius")
plt.title("Bootsrap for radius (Quatar-2b)")
plt.savefig("Bootsrap for radius (Quatar-2b)(n="+str(n)+")",dpi=256)
plt.show()
SDrad = fit.SD(Params[0])
print 'Error on radius = ' , SDrad

print 'Mass:'
plt.hist(Params[1],bins=15, color = 'blueviolet')
plt.ylabel("Count")
plt.xlabel("Mass /Solar mass")
plt.title("Bootsrap for mass (Quatar-2b)")
plt.savefig("Bootsrap for mass (Quatar-2b)(n="+str(n)+")",dpi=256)
plt.show()
SDmass = fit.SD(Params[1])
print 'Error on mass = ' , SDmass


x = np.arange(timeR[0],timeR[len(timeR)-1],1/1e4)
y = model(x, m['x'][0], m['x'][1], m['x'][2], m['x'][3], m['x'][4])
#y0 = model(time, x0[0], x0[1], x0[2], x0[3], x0[4])

plt.plot(x,y,color='green',lw=1.5)
plt.errorbar(timeR, relFlux, yerr=relFluxErr,linestyle='none',marker='x',capthick=0.3, mew=0.6,elinewidth=0.4, ms=4,mec="indigo",ecolor="blue") 

plt.ylabel("Relative flux")
plt.xlabel("Adjusted time /day")
plt.title("Eclipse from transit of Qatar-2b")
plt.savefig("Qatar-2b transit light curve with model curve",dpi=256)

plt.show()

if False:
    import datetime
    dt = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
    
    if not open('Qatar-2b value.txt', 'a'):
        open('Qatar-2b value.txt', 'w')
    
    file = open('Qatar-2b value.txt', 'a')
    file.write('Data from run completed on ' + dt +'\n')
    file.write('Planet radius = ' + str(m['x'][0])+'\n')
    file.write('Planet mass = ' + str(m['x'][1])+'\n')
    file.write('Planet t0 = ' + str(m['x'][2])+'\n')
    file.write('Scale = ' + str(m['x'][3])+'\n')
    file.write('Offset = ' + str(m['x'][4])+'\n')
    file.write('For n = ' + str(n)+'\n')
    file.write('Error on radius = ' + str(SDrad)+'\n')
    file.write('Error on mass = ' + str(SDmass)+'\n')
    file.write("\n-#-#-#-#-#-#-#\n\n")
    file.close()
    
print "\a"