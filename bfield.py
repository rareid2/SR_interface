import numpy as np
import datetime as dt
import os
from constants_settings import *
from convert_coords import convert2
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, get_yearmiliday
import PyGeopack as gp
import random
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

def getBline(pos,ray_datenum,stopalt):

    makedate = ray_datenum.strftime('%Y%m%d')
    Date = int(makedate)
    ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

    x = pos[0]/R_E
    y = pos[1]/R_E
    z = pos[2]/R_E

    # uses T96, need to confirm this is constistent w raytracer
    T = gp.TraceField(x,y,z,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM',alt=stopalt,dir=-1)
    return T

# --------------------------------------------------------------------
# OLD FUNCTION!! this requires runnign a ray to get the Bfield -- don't use this 
# use the next one for increased speed!!!!!!!!!!!!

# a helper func to run a ray to get the most accuarate Bfield 

# use for defining wavenormals

# inputs! 
# start pos (one position array)
# ray_datenum (one time)
# rayfile directory (one directory)

# if select_random = True
# thetas and phis are passed in as list of zeros of length nrays
# hemimult = 1 for north, -1 for south

# if select_random=False
# thetas and phis are passed in as list of length nrays where each pair is the desired polar and azimuthal angle to launch from

# output is converted_dirs (unit vector w directions) the resonance cone angle, and the thetas and phis passed in

def getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, hemimult, md, freq, select_random=False):
    positions = ray_start
    directions = [(0,0,0)]
    freqs = [freq]
    # going to run a ray real quick
    single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md, runmodeldump=False)

    # Load all the rayfiles in the output directory
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []
    for filename in file_titles:
        if '.ray' in filename and str(md) in filename and str('Main') in filename:
            print(filename)
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
    # this was used for conjunction analysis between DSX and VPM
    if select_random == True:
        nrays = len(thetas)
        thetas = []
        phis = []
        resangle_deg = resangle *180/np.pi

        for n in range(0,nrays):
            # sample theta as concentric circles around the z axis, max at resonance angle
            thetas.append((random.random()*(resangle_deg-1)))
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
            zsign = hemimult

            cone_vec = cone_vec/np.linalg.norm(cone_vec)
            converted_dirs.append(zsign*cone_vec)
            #ax.quiver(0,0,0,zsign*cone_vec[0],zsign*cone_vec[1],zsign*cone_vec[2],length=1)
        #ax.quiver(0,0,0,Bunit[0],Bunit[1],Bunit[2],length=1.25)
        #ax.axes.set_xlim3d(left=0, right=1.5) 
        #ax.axes.set_ylim3d(bottom=0, top=1.5) 
        #ax.axes.set_zlim3d(bottom=0, top=1.5) 
        #plt.show()
        #plt.close()

    # add theta and phi as desired
    else:
        for theta, phi in zip(thetas,phis):
            new_dir = [sph_dir[0][0],sph_dir[0][1]+theta, sph_dir[0][2]+phi]
            converted_dir = convert2([new_dir], ray_datenum, 'SM', 'sph', ['Re','deg','deg'], 'SM', 'car', ['Re','Re','Re']) 
            
            converted_dirs.append(converted_dir[0])

    return converted_dirs, resangle, thetas, phis # returns unit vector of directions corresponding to input theta and phi vals
# ---------------------------------------------------------------------------------
# this is the updated version of getting the magnetic field from the original src code
# WITHOUT running any rays - it is 100000 times faster, please use this one!

# inputs
# ray_datenum - datetime object when the rays are to be run
# positions - list of positions in SM cartesian in meters to sample Bfield at


def getBdir_src(ray_datenum, positions, thetas, phis):

    yearday, milliseconds_day = get_yearmiliday(ray_datenum)

    # had to add this in to fix issues with current directory
    cwd = os.getcwd()
    os.chdir(raytracer_dir)
    # make sure there is a model dumps folder
    ray_out_dir = cwd + '/modeldumps'

    model_outfile = os.path.join(ray_out_dir, 'bmodel_dump.dat')
    model_inpfile = os.path.join(ray_out_dir, 'bmodel_input.in')

    # Write input positions file
    f = open(model_inpfile, 'w')

    # Go through list of positions, write a new line for each
    for pos0 in positions:
        f.write('%1.15e %1.15e %1.15e\n' % (
            pos0[0], pos0[1], pos0[2]))
    f.close()

    # get number of rays for memory allocation
    nrays = len(positions)

    cmd = './getBfield ' + \
        '--filename="%s" ' % (model_outfile) + \
        '--infilename="%s" ' % (model_inpfile) + \
        ' --yearday=%s --milliseconds_day=%d ' % (yearday, milliseconds_day) + \
        '--use_igrf=%g --use_tsyganenko=%g ' % (use_IGRF, use_tsyg) + \
        '--tsyganenko_Pdyn=%g ' % (Pdyn) + \
        '--tsyganenko_Dst=%g ' % (Dst) + \
        '--tsyganenko_ByIMF=%g ' % (ByIMF) + \
        '--tsyganenko_BzIMF=%g ' % (BzIMF) + \
        '--tsyganenko_W1=%g ' % (W[0]) + \
        '--tsyganenko_W2=%g ' % (W[1]) + \
        '--tsyganenko_W3=%g ' % (W[2]) + \
        '--tsyganenko_W4=%g ' % (W[3]) + \
        '--tsyganenko_W5=%g ' % (W[4]) + \
        '--tsyganenko_W6=%g ' % (W[5]) + \
        '--nrays=%g' % nrays
    
    # execute the command
    os.system(cmd)

    # Move back to the working directory
    os.chdir(cwd)

    # save the outputs - converted directions
    converted_dirs = []

    # intermediate output
    sph_dirs = []   

    # open the file and save stuff
    f = open(model_outfile,'r')
    three_count = 0
    # read output file with all the directions
    for line in f.readlines():
        # comes out line by line, need to read in directions one at a time
        if three_count==0:
            B0 = []
        if three_count<2:
            B0.append(float(line))
            three_count+=1
        else:
            B0.append(float(line))
            three_count=0

            # reset, that means we have the direc
            Bmag = np.sqrt(B0[0]**2 + B0[1]**2 + B0[2]**2)
            Bunit = [B0[0]/Bmag, B0[1]/Bmag, B0[2]/Bmag]
            # now we have Bunit in SM car
            # let's put it in spherical (easier for changing the wavenormal)
            sph_dir = convert2([Bunit], ray_datenum, 'SM', 'car', ['Re','Re','Re'], 'SM', 'sph', ['Re','deg','deg'])     
            sph_dirs.append(sph_dir[0])
    f.close()

    # now change to new directions with thetas and phis
    for sph_dir, theta, phi in zip(sph_dirs,thetas,phis):
        new_dir = [sph_dir[0],sph_dir[1]+theta, sph_dir[2]+phi]
        # convert back to SM car
        converted_dir = convert2([new_dir], ray_datenum, 'SM', 'sph', ['Re','deg','deg'], 'SM', 'car', ['Re','Re','Re']) 
        
        converted_dirs.append(converted_dir[0])

    return converted_dirs