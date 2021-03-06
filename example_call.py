import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from satellites import sat
from bfield import getBdir
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, read_damp_simple
from ray_plots import plotray2D, plotrefractivesurface, plot_density_alongpath

# example call to ray tracer!
# --------------------------------------- set up ------------------------------------------
rayfile_directory = '/media/rileyannereid/DEMETER/SR_output' # store output here

# FIRST, navigate to constants_settings and make sure the settings are correct for the run

# let's look at a conjunction between DSX and VPM:
# use the datetime package to define the start time -- make sure to use UTC timezone
ray_datenum = dt.datetime(2014,1,1,12,0, tzinfo=dt.timezone.utc)

# we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # provide NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

# now we have ray start point in the correct coordinates (SM cartesian in m)
ray_start = dsx.pos

# Which plasmasphere model should we run?
#   1 - Legacy (Ngo) model
#   6 - Simplified GCPM from Austin Sousa's thesis
#   7 - New! Diffusive Equilibrium AT64ThCh (see docs)
md = 1

# how many rays? 
nrays = 1 # how many rays -- THIS MUST BE EQUAL IN LENGTH TO THETAS

# next, define the direction of the ray
# this step will actually run the raytracer to sample the Bfield correctly
# theta = 0 goes north, theta=180 goes south
thetas = [0 for i in range(nrays)] # go south
directions, ra, thetas, phis = getBdir(ray_start, ray_datenum, rayfile_directory, thetas, np.zeros(len(thetas)), md)

freq = 2.3e3 # Hz

positions = [ray_start[0] for n in range(nrays)]
freqs = [freq for n in range(nrays)]

# --------------------------------------- run ---------------------------------------------
# time to run is about 1 sec every 10 rays 
single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md, runmodeldump=True)

"""
# OR run parallel at different times -- use parallel_run_rays and input a list of times, and LIST OF LISTS with positions, 
# directions, and frequencies
"""

# ------------------------------------- output -------------------------------------------------
# that's it! let's look at output

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []
damplist = []
for filename in file_titles:
    if '.ray' in filename and str(md) in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=True)