import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from constants_settings import *
from bfield import Bfieldinfo, trace_fieldline_ODE

#from spc_coordinates import create_spc, convert_spc
from libxformd import xflib
xf = xflib.xflib(lib_path='libxformd/libxformd.so')

# -------------------------------- TX CLASS -----------------------------------------
# just a nice way to mangage the transmitter launches

# pos input must be spherical GEO coords

# see end of script for example call
# -------------------------------------------------------------------------------------

class vlf_tx:
    def __init__(self, time=None, pos=None, freq=None):
        self.time = None
        self.pos = None
        self.freq = None

    def tracepos_up_fieldline(self, start_alt): # alt in meters
        if self.time == None:
            print('need time first')
            return
        if self.pos == None:
            print('need pos first')
            return

        if float(check_hemi.lati) > 0:
            hemis = 'south' # go opposite the tx
        else:
            hemis = 'north'

        tx_b = Bfieldinfo()
        tx_b.time = self.time 
        tx_b.pos = self.pos

        trace_fieldline_ODE(tx_b, hemis, 'GEO', 'sph', ['m','m','m'])

        bline = tx_b.fieldline
        trace_alt = R_E + start_alt

        for pos_crs in bline:
            if float(pos_crs.radi) > trace_alt:
                traced_pos = pos_crs
                break

        tx_crs_traced = convert_spc(traced_pos, self.time,'SM','car',['m','m','m'])

        return tx_crs_traced

# ---------------------------------------------------------------------

# Example Call! 
# let's look at a TX at a specific date and time
#ray_datenum = dt.datetime(2020, 9, 14, 22, 55, tzinfo=dt.timezone.utc)

# first, where is the transmitter on earth -- MUST be in GEO sph coordinates in m (or Re), deg, deg
#tx_loc = [1, 33.2385, -106.3464]

#wsmr = vlf_tx()
#wsmr.time = ray_datenum
#wsmr.pos = tx_loc
#wsmr.freq = 14.1e3

# returns in SM coords
#tx_crs_traced = wsmr.tracepos_up_fieldline(start_alt=1000e3)

# ---------------------------------------------------------------------