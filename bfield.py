import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
import PyGeopack as gp

# uses PyGeopack (same as raytracer but its still a lil diff)
# pos MUST be in SM car in m
# returns line in SM car in m

# output T object:
"""
x 	x coordinate along the field trace(s)
y 	y coordinate along the field trace(s)
z 	z coordinate along the field trace(s)
Bx 	x component of the magnetic field along the trace(s)
By 	y component of the magnetic field along the trace(s)
Bz 	z component of the magnetic field along the trace(s)
nstep 	number of steps along the trace(s)
GlatN 	Geographic latitude of the northern footprint(s)
GlatS 	Geographic latitude of the southern footprint(s)
MlatN 	Magnetic latitude of the northern footprint(s)
MlatS 	Magnetic latitude of the southern footprint(s)
GlonN 	Geographic longitude of the northern footprint(s)
GlonS 	Geographic longitude of the southern footprint(s)
MlonN 	Magnetic longitude of the northern footprint(s)
MlonS 	Magnetic longitude of the southern footprint(s)
GltN 	Geographic local time of the northern footprint(s)
GltS 	Geographic local time of the southern footprint(s)
MltN 	Magnetic local time of the northern footprint(s)
MltS 	Magnetic local time of the southern footprint(s)
Lshell 	L-shell of the field line(s) at the equator
MltE 	Magnetic local time of the equatorial footprint(s)
FlLen 	Field line length in planetary radii
R 	R = sqrt(x**2 + y**2 + z**2)
"""

def getBline(pos, ray_datenum, stopalt):

    makedate = ray_datenum.strftime('%Y%m%d')
    Date = int(makedate)
    ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

    x = pos[0]/R_E
    y = pos[1]/R_E
    z = pos[2]/R_E

    # uses T96, need to confirm this is constistent w raytracer
    T = gp.TraceField(x,y,z,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM')
    return T

# --------------------------------------------------------------------
# a helper func to run a ray to get the most accuarate Bfield 
# use for defining wavenormals

def getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis):
    positions = ray_start
    directions = [(0,0,0)]
    freqs = [5e3]

    single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory)

    # Load all the rayfiles in the output directory
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []
    for filename in file_titles:
        if '.ray' in filename:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))

    for r in raylist:
        B0 = [r['B0'].x[0], r['B0'].y[0], r['B0'].z[0]]
        # vec in T in SM car coordinates
        # create unit vector 
        Bmag = np.sqrt(r['B0'].x[0]**2 + r['B0'].y[0]**2 + r['B0'].z[0]**2)
        Bunit = [r['B0'].x[0]/Bmag, r['B0'].y[0]/Bmag, r['B0'].z[0]/Bmag]

    # now we have Bunit in SM car
    # let's put it in spherical (easier for changing the wavenormal)
    sph_dir = convert2([Bunit], ray_datenum, 'SM', 'car', ['Re','Re','Re'], 'SM', 'sph', ['Re','deg','deg'])     
    
    converted_dirs = []
    # add theta and phi
    for theta, phi in zip(thetas,phis):
        new_dir = [sph_dir[0][0],sph_dir[0][1]+theta, sph_dir[0][2]+phi]
        converted_dir = convert2([new_dir], ray_datenum, 'SM', 'sph', ['Re','deg','deg'], 'SM', 'car', ['Re','Re','Re']) 
        converted_dirs.append(converted_dir[0])
    return converted_dirs # returns unit vector of directions corresponding to input theta and phi vals