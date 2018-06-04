#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 18:46:56 2018

@author: twsee
"""

import xarray
import glob
import pandas as pd
import numpy as np
from wrf import getvar
def gatherfiles(prefix):
    return glob.glob(prefix+'*.nc')
    
files = gatherfiles('wrfout')
ncfile = xarray.Dataset("wrfout_d01_2016-08-02_12:00:00.nc")

#power =  Cp 1/2 rho A V**3
#power is power output in kilowats
#Cp = Maximum power coefficient, ranging from 0.25 to 0.45, dimension less (theoretical maximum = 0.59)
#ρ = Air density, lb/ft3
#A = Rotor swept area, ft2 or π D2/4 (D is the rotor diameter in ft, π = 3.1416) rotor diameter seems to be around 125m according to a site, will make a range of it. 
#V = Wind speed, m/s

#to get height(altitude) i need to calculate it by (PB + PHB) / g
r = 52.0                            #meters, rotorlength
Area = np.pi * r**2
cp = 0.4                            #unitless       
density = 1.23                      #kg/m^3
press = getvar(ncfile, 'pressure')  #hPa
T = getvar(ncfile,'T')              #K
ua = getvar(ncfile, 'ua')           #m/s
va = getvar(ncfile, 'va')           #m/s
pressground, uaground, vaground,Tground = press[0,:,:],ua[0,:,:], va[0,:,:], T[0,:,:]
#rho = pressground/(R*T)
totalwind = np.sqrt(uaground**2 + vaground**2)

energyproduction = (0.5 * Area * totalwind**3 * cp)/10**6 #megawatts output


def gathervariables(dataframe):
    dataframe = xarray.open_dataset(dataframe)
    uwind = dataframe['U']
    vwind = dataframe['V']
    U = uwind.to_dataframe()
    V = vwind.to_dataframe()
    return U,V