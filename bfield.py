import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
#import PyGeopack as gp
from random import randrange, uniform
import matplotlib.pyplot as plt

# uses PyGeopack (same as raytracer but its still a lil diff)
# pos MUST be in SM car in m
# returns line in SM car in Re

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

def getBline(pos, ray_datenum,stopalt):

    makedate = ray_datenum.strftime('%Y%m%d')
    Date = int(makedate)
    ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

    x = pos[0]/R_E
    y = pos[1]/R_E
    z = pos[2]/R_E

    # uses T96, need to confirm this is constistent w raytracer
    T = gp.TraceField(x,y,z,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM',alt=stopalt)
    return T

# --------------------------------------------------------------------
# a helper func to run a ray to get the most accuarate Bfield 
# use for defining wavenormals

def getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, md, select_random=False):
    positions = ray_start
    directions = [(0,0,0)]
    freqs = [15e3]
    # going to run a ray real quick
    single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)

    # Load all the rayfiles in the output directory
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []
    for filename in file_titles:
        if '.ray' in filename and str(md) in filename:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))

    # get b direction for this ray
    for r in raylist:
        B0 = [r['B0'].x[0], r['B0'].y[0], r['B0'].z[0]]
        # vec in T in SM car coordinates
        # create unit vector 
        Bmag = np.sqrt(r['B0'].x[0]**2 + r['B0'].y[0]**2 + r['B0'].z[0]**2)
        Bunit = [r['B0'].x[0]/Bmag, r['B0'].y[0]/Bmag, r['B0'].z[0]/Bmag]

    # now we have Bunit in SM car
    # let's put it in spherical (easier for changing the wavenormal)
    sph_dir = convert2([Bunit], ray_datenum, 'SM', 'car', ['Re','Re','Re'], 'SM', 'sph', ['Re','deg','deg'])     

    # also return resonance angle, can be useful for initializing rays
    from ray_plots import stix_parameters
    R, L, P, S, D = stix_parameters(r, 0, r['w']) # get stix params for initial time point
    resangle = np.arctan(np.sqrt(-P/S))
    
    converted_dirs = []
    # if select random was chosen, thetas and phis are passed in as list of zeros of length nrays
    if select_random == True:
        nrays = len(thetas)
        hemi_mult = thetas[0]
        thetas = []
        phis = []
        resangle_deg = resangle *180/np.pi

        for n in range(0,nrays):
            # sample theta as concentric circles around the z axis, max at resonance angle
            thetas.append((random.random()*(resangle_deg-3)))
            # uniform azimuth around the z axis
            phis.append(random.random()*360)

        if Bunit[0] == 0 or Bunit[2] == 0:
            r1 = [1,(-1*Bunit[0]-Bunit[2])/Bunit[1],1]
        else:
            r1 = [1,1,(-1*Bunit[1]-Bunit[0])/Bunit[2]]

        r1 = np.array(r1)/np.linalg.norm(np.array(r1))
        r2  = np.cross(r1,Bunit)
        T_rotate = np.column_stack((r1,r2,Bunit))

        #ax = plt.axes(projection='3d')
        for th,ph in zip(thetas,phis):
            r = 1/(np.cos(th*D2R))
            cone_vec = np.array([r*np.sin(th*D2R)*np.cos(ph*D2R),r*np.sin(th*D2R)*np.sin(ph*D2R),r*np.cos(th*D2R)]) 
            cone_vec = np.matmul(T_rotate,np.transpose(cone_vec))
            if hemi_mult == 180:
                zsign = -1
            else:
                zsign = 1

            cone_vec = cone_vec/np.linalg.norm(cone_vec)
            converted_dirs.append(zsign*cone_vec)
            #ax.plot3D([0,cone_vec[0]], [0,cone_vec[1]], [0,zsign*cone_vec[2]])
        #ax.plot3D([0,Bunit[0]],[0,Bunit[1]],[0,Bunit[2]],'k')
        #plt.show()
        #plt.close()

    # add theta and phi as desired
    else:
        for theta, phi in zip(thetas,phis):
            new_dir = [sph_dir[0][0],sph_dir[0][1]+theta, sph_dir[0][2]+phi]
            converted_dir = convert2([new_dir], ray_datenum, 'SM', 'sph', ['Re','deg','deg'], 'SM', 'car', ['Re','Re','Re']) 
            
            converted_dirs.append(converted_dir[0])

    return converted_dirs, resangle, thetas, phis # returns unit vector of directions corresponding to input theta and phi vals