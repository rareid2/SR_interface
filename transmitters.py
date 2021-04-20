import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from constants_settings import *
from bfield import getBline
from convert_coords import convert2
from libxformd import xflib
xf = xflib.xflib(lib_path='libxformd/libxformd.so')

# -------------------------------- TX CLASS -----------------------------------------
# just a nice way to mangage the transmitter launches

# see end of script for example call
# -------------------------------------------------------------------------------------

class vlf_tx:
    def __init__(self, time=None, pos=None, freq=None):
        self.time = None
        self.pos = None
        self.freq = None

    def tracepos_up_fieldline(self, start_alt, crs, carsph, units): # alt in meters
        if self.time == None:
            print('need time first')
            return
        if self.pos == None:
            print('need pos first')
            return

        # bfield needs input in SM car, m
        bpos = convert2([self.pos], [self.time], crs, carsph, units, 'SM','car', ['m','m','m'])

        T = getBline(bpos[0],self.time,100) # returned in SM coords car Re
        
        # repack and convert
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        LT = [self.time for i in range(len(T_repack))]
        
        tx_crs_traced = convert2(T_repack, LT, 'SM','car',['Re','Re','Re'], crs, carsph, units)

        # to find location at start alt, need to interpolate
        start_alt = ((start_alt*1e3) / R_E) + 1 # now in RE

        # why did this throw an error? 
        if carsph != 'sph': #need sph for this to work, add more handling later
            print('error!')
        else:
            for ti,tx in enumerate(tx_crs_traced):
                if start_alt - tx[0] < 0.01:
                    p2 = tx[0]
                    savelong = tx[2]
                    save_ind = ti
                    break
            p1 = tx_crs_traced[ti-1][0]
            
            from scipy.interpolate import interp1d  
            xs = [tx[0] for tx in tx_crs_traced[ti-2:ti+2]]
            ys = [tx[1] for tx in tx_crs_traced[ti-2:ti+2]]
            f2 = interp1d(xs, ys, kind='cubic') # interpolate over area we are looking for

            # new x to interpolate over
            new_xs = np.linspace(p1,p2,num=100,endpoint=True)
            interpcoords = f2(new_xs) # interpolated latitudes
            
            for xi, x in enumerate(new_xs):
                if np.abs(x-start_alt) < 0.0001: # find closest pt
                    savenx = x
                    saveny = interpcoords[xi]
            
            # assume longitude does not change
            traced_pos = [savenx, saveny, savelong]

        return tx_crs_traced, traced_pos

# ---------------------------------------------------------------------

# Example Call! 
# let's look at a TX at a specific date and time
#ray_datenum = dt.datetime(2020, 9, 14, 22, 55, tzinfo=dt.timezone.utc)

# first, where is the transmitter on earth -- set crs carsph units to be input and output
#crs='GEO'
#carsph='sph'
#units=['Re','deg','deg'] # NEEDS TO BE RE
#tx_loc = [1, 33.2385, -106.3464]

#wsmr = vlf_tx()
#wsmr.time = ray_datenum
#wsmr.pos = tx_loc
#wsmr.freq = 14.1e3

# find position traced up fieldline
#tx_crs_traced, traced_pos = wsmr.tracepos_up_fieldline(1000,crs,carsph,units) # input start alt in km
# ---------------------------------------------------------------------