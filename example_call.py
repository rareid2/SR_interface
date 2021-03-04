import numpy as np
import datetime as dt
import os

from constants_settings import *
from satellites import sat
from bfield import getBdir
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
from ray_plots import plotray2D, plotrefractivesurface

# example call to ray tracer!
# FIRST, navigate to constants_settings and make sure the settings are correct for the run

# let's look at a conjunction between DSX and VPM:
# use the datetime package to define the start time -- make sure to use UTC timezone
ray_datenum = dt.datetime(2020,8,17,21,22,tzinfo=dt.timezone.utc)

# we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # provide NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

# now we have ray start point in the correct coordinates (SM cartesian in m)
ray_start = dsx.pos

vpm = sat()             # define a satellite object
vpm.catnmbr = 45120     # provide NORAD ID
vpm.time = ray_datenum  # set time
vpm.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag
vpm.propagatefromTLE(sec=0, orbit_dir='future', crs='GEO', carsph='sph', units=['m','deg','deg'])
print(vpm.pos)

# next, define the direction of the ray
# this step will actually run the raytracer to sample the Bfield correctly
# returns the Bfield unit vector at the start point in SM car

# thetas and phis are vectors corresponding to starting directions
# theta = increases the polar angle
# phi = increases the azimuth angle
rayfile_directory = '/home/rileyannereid/workspace/SR_output'

# returns a vector of directions (thetas and phis must be same length) 
directions = getBdir(ray_start, ray_datenum, rayfile_directory, thetas=[0], phis=[0])

# run at a single time -- use run_rays and input a list of positions, directions, and freqs (ALL SAME LENGTH)
# generates one input file and one output file with all rays in it

nrays = 1 # how many rays -- THIS MUST BE EQUAL IN LENGTH TO THETAS AND PHIS
freq = 8.83e3  # Hz

# simplest call
positions = [ray_start[0] for n in range(nrays)]
freqs = [freq for n in range(nrays)]

# time to run is about 1 sec every 10 rays 
single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory)

"""
# OR run parallel at different times -- use parallel_run_rays and input a list of times, and LIST OF LISTS with positions, 
# directions, and frequencies

# let's say we want to re-run this every 10 seconds in time for a minute
tvec = [ray_datenum + dt.timedelta(seconds=i) for i in range(0,60,10)]
positions_list = [positions for i in range(len(tvec))]
directions_list = [directions for i in range(len(tvec))]
freqs_list = [freqs for i in range(len(tvec))]
directory_list = [rayfile_directory for i in range(len(tvec))]

parallel_run_rays(tvec, positions_list, directions_list, freqs_list, directory_list)
"""

# -------------------------------------------------------------------------
# that's it! let's look at output

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []
for filename in file_titles:
    if '.ray' in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', units=['Re','Re','Re'])
plotrefractivesurface(ray_datenum, raylist[0], ray_out_dir)