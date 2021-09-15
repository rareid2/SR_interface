import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from satellites import sat
from bfield import getBdir
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, read_damp_simple
from ray_plots import plotray2D, plotrefractivesurface, plot_plasmasphere_2D

# example call to ray tracer!
# --------------------------------------- set up ------------------------------------------
# STEP 1 - set the output directory
rayfile_directory = '/media/rileyannereid/DEMETER/SR_output' # store output here

# STEP 2 - navigate to constants_settings.py and make sure the settings are correct for the run

# STEP 3 - set the date and time
# use the datetime package to define the start time -- make sure to use UTC timezone
ray_datenum = dt.datetime(2020,6,6,19,55,0, tzinfo=dt.timezone.utc)

# STEP 4 - set the start position
# here's how to get it from satellites: 
# we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # provide NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

# now we have ray start point in the correct coordinates (SM cartesian in m)
ray_start = dsx.pos

# here's how to get it from a start location in GEO coords 
start_lat = 50
start_lon = 0
start_alt = 1000e3 # m
ray_start = convert2([[start_alt+R_E,start_lat,start_lon]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])

# STEP 5 - set freq in Hz
freqs = [8e3] # Hz

# STEP 6 - which plasmasphere model should we run?
#   6 - Simplified GCPM from Austin Sousa's thesis
#   7 - New! Diffusive Equilibrium AT64ThCh (see docs)
md = 7

# STEP 7 - how many rays? 
nrays = 1

# STEP 8 - in what direction should the ray point? 
# set theta as the initial wavenormal angle from Bdir
# range from -90 to 90 deg will go north, 90 to 270 will go south
thetas = [180]
# phi is in the longitudanal direction
phis = [0]
# hemimult = 1 for noth, -1 for south
hemimult = -1
directions, ra, thetas, phis = getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, hemimult, md, freqs[0])

# STEP 9 - run, time to run is about 1 sec every 10 rays 
single_run_rays(ray_datenum, ray_start, directions, freqs, rayfile_directory, md, runmodeldump=False)

# ------------------------------------- output -------------------------------------------------
# STEP 10 - that's it! let's look at output

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files
raylist = []
for filename in file_titles:
    if '.ray' in filename and str(md) in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

# example plot
plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', ['Re','Re','Re'], md, show_plot=True)