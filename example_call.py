import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from satellites import sat
from bfield import getBdir, getBdir_src
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, read_damp_simple, get_yearmiliday
from ray_plots import plotray2D, plotrefractivesurface, plot_plasmasphere_2D

# example call to ray tracer!
# --------------------------------------- set up ------------------------------------------
# STEP 1 - set the output directory
rayfile_directory = '/home/rileyannereid/workspace/2020-06-01_12_00_00_NWCtest' # store output here

# STEP 2 - navigate to constants_settings.py and make sure the settings are correct for the run

# STEP 3 - set the date and time
# use the datetime package to define the start time -- make sure to use UTC timezone
ray_datenum = dt.datetime(2020,6,1,12,0,0, tzinfo=dt.timezone.utc)

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
#start_lat = -21.82
#start_lon =  114.17
#start_alt = 0 # m
#ray_start = convert2([[start_alt+R_E,start_lat,start_lon]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])

# STEP 5 - set freq in Hz
freqs = [19.88e3] # Hz

# STEP 6 - which plasmasphere model should we run?
#   6 - Simplified GCPM from Austin Sousa's thesis
#   7 - New! Diffusive Equilibrium AT64ThCh (see docs)
md = 7

# STEP 7 - how many rays? 
nrays = 1

# STEP 8 - in what direction should the ray point? -- there are 3 methods, they all seems to be correct within a few deg
# set theta as the initial wavenormal angle from Bdir - this is a cone around the Bdirection with the angle theta
# range from -90 to 90 deg will go north, 90 to 270 will go south
thetas = [0]
# phi is the azimuth (on the surface of the cone -- see Stanford Raytracer gihtub for a visualization)
phis = [0]

directions = getBdir_src(ray_datenum, ray_start, thetas, phis)

# STEP 9 - run, time to run is about 1 sec every 10 rays 
single_run_rays(ray_datenum, ray_start, directions, freqs, rayfile_directory, md, runmodeldump=False)

# ------------------------------------- output -------------------------------------------------
# STEP 10 - that's it! let's look at output

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
file_titles = os.listdir(ray_out_dir)
file_titles.sort()
# create empty lists to fill with ray files
raylist = []
for filename in file_titles:
    if '.ray' in filename and str(md):
        print(filename)
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

# example plot
plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', ['Re','Re','Re'], md, show_plot=True, plot_wna=True)


""" # heres a loop to 'trace' a field line with a given step size
stopalt = 1000 # km
dstep = 10 # km
n_steps = int(stopalt / dstep)
for i in range(n_steps):
    directions = getBdir_src(ray_datenum, ray_start_surface, thetas, phis)

    newdir = dstep*1e3 * directions[0] # in meters 
    ray_start = [ray_start_surface[0][0] + newdir[0], ray_start_surface[0][1] + newdir[1], ray_start_surface[0][2] + newdir[2]]
    # update
    ray_start_surface = [ray_start]
"""