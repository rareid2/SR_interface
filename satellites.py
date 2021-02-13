import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import requests
import configparser
import json
import time
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv

from constants_settings import *
from convert_coords import convert2
from libxformd import xflib
xf = xflib.xflib(lib_path='libxformd/libxformd.so')

from bfield_for import bfieldinfo

# -------------------------------- SAT CLASS -----------------------------------------
# just a nice way to mangage the satellites (VPM, DSX...)

# provide catnmbr (44344 or 45120 for DSX and VPM respectively)
# provide time in UTC

# return TLE, orbit (provide desired length to propagate the orbit and direction of propagation)
# and get velocity interesting for doppler

# see end of script for example call
# -------------------------------------------------------------------------------------

class sat:
    def __init__(self, TLE=None, time=None, pos=None, vel=None):
        self.catnmbr = None # catnmbr is NOARD ID
        self.TLE = None
        self.time = None
        self.pos = None
        self.vel = None

    def getTLE_ephem(self):
        if self.catnmbr == None:
            print('need catnmbr first')
            return
        if self.time == None:
            print('need time first')
            return

        class MyError(Exception):
            def __init___(self,args):
                Exception.__init__(self,"my exception was raised with arguments {0}".format(args))
                self.args = args

        uriBase                = "https://www.space-track.org"
        requestLogin           = "/ajaxauth/login"
        requestCmdAction       = "/basicspacedata/query" 
        requestFindSat         = "/class/tle/NORAD_CAT_ID/"+str(self.catnmbr)+"/"

        # Use configparser package to pull in the ini file (pip install configparser)
        config = configparser.ConfigParser()
        config.read("./SLTrack.ini")
        configUsr = config.get("configuration","username")
        configPwd = config.get("configuration","password")
        configOut = config.get("configuration","output")
        siteCred = {'identity': configUsr, 'password': configPwd}

        # use requests package to drive the RESTful session with space-track.org
        with requests.Session() as session:
            # run the session in a with block to force session to close if we exit

            # need to log in first. note that we get a 200 to say the web site got the data, not that we are logged in
            resp = session.post(uriBase + requestLogin, data = siteCred)
            if resp.status_code != 200:
                raise MyError(resp, "POST fail on login")

            # this query picks up all Starlink satellites from the catalog. Note - a 401 failure shows you have bad credentials 
            resp = session.get(uriBase + requestCmdAction + requestFindSat)
            if resp.status_code != 200:
                raise MyError(resp, "GET fail on request for satellite")

            # use the json package to break the json formatted response text into a Python structure (a list of dictionaries)
            retData = json.loads(resp.text)
            epoch_dt = [] # save epoch tiemdelta
            for n in range(len(retData)):
                tle_time = dt.datetime.strptime(retData[n]['EPOCH'],'%Y-%m-%d %H:%M:%S').replace(tzinfo=dt.timezone.utc)
                tle_us = dt.timedelta(microseconds=float(retData[n]['EPOCH_MICROSECONDS']))
                tle_fulltime = tle_time + tle_us
                epoch_dt.append(np.abs((self.time - tle_fulltime).total_seconds()))
            
            epoch_n = epoch_dt.index(min(epoch_dt))
            # get the data!
            self.TLE = [retData[epoch_n]['TLE_LINE1'], retData[n]['TLE_LINE2']]

        session.close()
        print("retrieved TLE from SpaceTrak for " + retData[epoch_n]['EPOCH'])
    
    def propagatefromTLE(self, sec, orbit_dir, crs, carsph, units):
        if self.TLE == None:
            print('need TLE')
            return

        if orbit_dir == 'future':
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(0,sec+1)]
        elif orbit_dir == 'past':
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(-sec-1,0)]
        else:
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(-sec-1,sec+1)] 

        # store pos and vel data
        pos_array = np.zeros((len(dt_array),3))
        vel_array = np.zeros((len(dt_array),3))

        for ti, tt in enumerate(dt_array):
            satellite = twoline2rv(self.TLE[0], self.TLE[1], wgs84)
            position, velocity = satellite.propagate(tt.year,tt.month,tt.day,tt.hour,tt.minute,tt.second)
            pos_array[ti,:] = position
            vel_array[ti,:] = velocity

        # finally, convert (starts from GEI coords in km cartesian)
        converted_coords = convert2(pos_array, dt_array, 'GEI', 'car', ['km','km','km'], crs, carsph, units)
        self.pos = converted_coords
        self.time = dt_array # update time as well 

    def getB(self, hemis, currentcrs, currentcarsph, currentunits):
        if currentcrs != 'GEO' or currentcarsph != 'sph' or currentunits != ['km','deg','deg']:
            newcoords = convert2(self.pos, self.time, crs, carsph, units, 'GEO', 'sph', ['km','deg','deg'])
        else:
            newcoords = self.pos
        
        # call bfield class and add pos and time (now in GEO sph coords)
        b = bfieldinfo()
        b.pos = newcoords
        b.time = self.time
        b.getBfield()

        # defined desired direction 
        if hemis == 'south':
            dir = -1
        else:
            dir = 1

        self.Bdir = b.Bunit * dir



# ---------------------------------------------------------------------

# Example Call! 
# define an obj called dsx
dsx = sat()

# set the NORAD ID for dsx
dsx.catnmbr = 44344

# set the time we want orbit info
dsx.time = dt.datetime(2020,5,6,11,10, tzinfo=dt.timezone.utc)

# get TLE
dsx.getTLE_ephem()

# propagate for 10 seconds
psec = 10

# propagate into the future
# alternatively, use 'past' for previous
# or 'both' for propagating in both directions in time
pdir = 'past'

# select desired coordinate system only GEO and SM supported rn
crs = 'GEO'
carsph = 'sph'
units = ['Re','deg','deg']

# propagate
dsx.propagatefromTLE(int(psec), pdir, crs, carsph, units)
dsx.getB('north', crs, carsph, units)
print(dsx.Bdir)
# look @ updated position vec