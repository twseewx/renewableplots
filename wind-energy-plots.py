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
from wrf import to_np, getvar, smooth2d, get_basemap, latlon_coords
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap


def gatherfiles(prefix):
    return glob.glob(prefix+'*.nc')
    
files = gatherfiles('wrfout')
ncfile = Dataset("wrfout_d01_2016-08-02_12:00:00.nc")

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
def gettimeanddates(file):
    split = file.split('_')
    date = split[2]
    hour = split[3].split(':')[0]
    return date, hour
def energyproduction(files, level):
    """
    Creates hourly and daily energy production outputs for wind a turbine
    output plots

    Parameters
    ----------
    arg1: list
        A 1D list of the file outputs from WRF
    arg2: integer
        level of the WRF module to run the evaluation on, typicall 0-2
        for low levels of the atmosphere. 

    Returns
    -------
    None
    """
    count = 0
    dailyenergy = 0
    for file in files:
        date, hour = gettimeanddates(file)
        ncfile = Dataset(file)
        cbarticks=np.arange(0.0,10.0,0.5)
        r = 52.0                            #meters, rotorlength
        Area = np.pi * r**2
        cp = 0.4                            #unitless       
        density = 1.23                      #kg/m^3
        press = getvar(ncfile, 'pressure')  #hPa
        T = getvar(ncfile,'T')              #K
        ua = getvar(ncfile, 'ua')           #m/s
        va = getvar(ncfile, 'va')           #m/s
        pressground, uaground, vaground,Tground = press[level,:,:],ua[level,:,:],va[level,:,:],T[level,:,:]
        bm = get_basemap(pressground)
        fig = plt.figure(figsize=(12,9))
        totalwind = np.sqrt(uaground**2 + vaground**2)
        energyproduction = (0.5 * density * Area * totalwind**3 * cp)/10**6#megawatts output
        lat,lon = latlon_coords(pressground)
        x,y = bm(to_np(lon),to_np(lat))
        bm.drawcoastlines(linewidth=0.25)
        bm.drawstates(linewidth=0.25)
        bm.drawcountries(linewidth=0.25)
        bm.contour(x, y, to_np(energyproduction), cbarticks, colors="black",vmin=0,vmax=10.0)
        bm.contourf(x, y, to_np(energyproduction), cbarticks,cmap = get_cmap('jet'),vmin=0,vmax=10.0)
        plt.colorbar(shrink=.62,ticks=cbarticks)
        plt.title('Hourly Energy Production for '+ date + ' ' + hour)
        plt.show()
        if count ==23:
            bm = get_basemap(pressground)
            fig = plt.figure(figsize=(12,9))
            totalwind = np.sqrt(uaground**2 + vaground**2)
            lat,lon = latlon_coords(pressground)
            x,y = bm(to_np(lon),to_np(lat))
            bm.drawcoastlines(linewidth=0.25)
            bm.drawstates(linewidth=0.25)
            bm.drawcountries(linewidth=0.25)
            dailyenergy = dailyenergy/24.0
            bm.contour(x, y, to_np(dailyenergy), cbarticks, colors="black",vmin=0,vmax=10.0)
            bm.contourf(x, y, to_np(dailyenergy), cbarticks,cmap = get_cmap('jet'),vmin=0,vmax=10.0)
            plt.colorbar(shrink=.62,ticks=cbarticks)
            plt.title('Dialy Average ' + yesterdaysdate)
            plt.show()
            dailyenergy=0
            count=0
        else:
            dailyenergy = energyproduction + dailyenergy
            print(dailyenergy)
            count = count+1
        yesterdaysdate=date
    #if count ==23:
    #    dailyenergy =0
    #dailyenergy = energyproduction + dailyenergy




def gathervariables(dataframe):
    dataframe = xarray.open_dataset(dataframe)
    uwind = dataframe['U']
    vwind = dataframe['V']
    U = uwind.to_dataframe()
    V = vwind.to_dataframe()
    return U,V