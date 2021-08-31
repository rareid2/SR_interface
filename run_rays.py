# import needed packages
import numpy as np
import os
import datetime as dt
from multiprocessing import Pool, cpu_count
import multiprocessing
from constants_settings import *
from run_model_dump import modeldump
from raytracer_utils import get_yearmiliday

def single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, mode, runmodeldump=False):

    # this is just a silly way to get the names of the workers
    fp = str(multiprocessing.current_process())[19:36]
    fp = fp.replace("'",'')
    if len(fp) < 17:
        fp = fp[:-1] + '0' + fp[-1]
    # make numbers double digit

    yearday, milliseconds_day = get_yearmiliday(ray_datenum)

    project_root = os.getcwd()  # grabs current full path

    # Create directory for inputs/outputs if doesn't already exist
    ray_out_dir = rayfile_directory + '/'+ dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S') #+ '/' + str(int(freqs[0]/10))

    if not os.path.exists(ray_out_dir):
        os.makedirs(ray_out_dir)

    # Set input file path
    ray_inpfile = os.path.join(ray_out_dir, 'ray_inpfile'+fp+'.txt')

    if runmodeldump:
        modeldump(ray_datenum, mode)

    # Set config file for mode 3
    mode3_interpfile = os.path.join(project_root, 'precomputed_grids',
                                                'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
    # Set config file for mode 4
    mode4_modelfile = os.path.join(project_root,
                                'precomputed_grids', 'precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

    # Write the ray input file
    f = open(ray_inpfile, 'w')

    # Go through list of positions, write a new ray for every direction and freq at each pos
    for pos0, dir0, fr in zip(positions, directions, freqs):
        w0 = fr * 2.0 * np.pi
        f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n' % (
            pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
    f.close()

    # GCPM model and damping code needs to be run in the same directory
    # as the binary file (and all the misc data files)
    cwd = os.getcwd()
    os.chdir(raytracer_dir)

    # Set output file path
    ray_outfile = os.path.join(ray_out_dir, 'ray_out_mode%d_%s.ray' % (mode,fp))
    damp_outfile = os.path.join(ray_out_dir, 'ray_out_mode%d.damp' % mode)

    # The base command -- with parameters common for all modes
    # what is outputper? 
    base_cmd = './raytracer --outputper=%d --dt0=%g --dtmax=%g' % (1, dt0, dtmax) + \
                ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g' % (t_max, root, fixedstep, maxerr) + \
                ' --modelnum=%d --maxsteps=%d --minalt=%d --inputraysfile="%s"' % (
                    mode, maxsteps, minalt, ray_inpfile) + \
                ' --outputfile="%s" --yearday=%s --milliseconds_day=%d' % (ray_outfile, yearday, milliseconds_day) + \
                ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g' % (use_tsyg, use_IGRF, Pdyn) + \
                ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g' % (Dst, ByIMF, BzIMF) + \
                ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g' % (W[0], W[1], W[2]) + \
                ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g' % (W[3], W[4], W[5])

    
    base_damp_cmd = './damping --inp_file "%s" --out_file "%s" ' % (ray_outfile, damp_outfile) + \
                    ' --Kp %g --AE %g' % (Kp, AE) + \
                    ' --yearday %s --msec %d' % (yearday, milliseconds_day) + \
                    ' --geom_factor=%d' % include_geom_factor
    
    if mode == 1:
        # Test the Ngo model
        damp_mode = 0

        ray_cmd = base_cmd + ' --ngo_configfile="%s"' % (configfile)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

    if mode == 2:
        # Test the full GCPM model
        damp_mode = 1

        ray_cmd = base_cmd + ' --gcpm_kp=%g' % (Kp)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

    if mode == 3:
        # Test the uniformly-sampled GCPM model
        damp_mode = 1

        ray_cmd = base_cmd + ' --interp_interpfile="%s"' % (mode3_interpfile)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

    if mode == 4:
        # Test the randomly-sampled GCPM model
        damp_mode = 1

        ray_cmd = base_cmd + ' --interp_interpfile=%s' % (mode4_modelfile) + \
                    ' --scattered_interp_window_scale=%d' % (scattered_interp_window_scale) + \
                    ' --scattered_interp_order=%d' % (scattered_interp_order) + \
                    ' --scattered_interp_exact=%d' % (scattered_interp_exact) + \
                    ' --scattered_interp_local_window_scale=%d' % (scattered_interp_local_window_scale)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

    if mode == 6:
        # Test the Simplified GCPM model
        MLT = 0
        fixed_MLT = 0 # Force the raytracer to stay in the meridonal plane?
        damp_mode = 0

        ray_cmd = base_cmd + ' --MLT="%g" --fixed_MLT=%g --kp=%g' % (MLT, fixed_MLT, Kp)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode
    
    if mode == 7:
        # Test the difusive equilibrium model from ARFL
        damp_mode = 1

        ray_cmd = base_cmd + ' --gcpm_kp=%g' % (Kp)
        damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

    # Run it!

    #print("------- Running mode %d -------" % mode)
    #print("Command is:")
    #print(ray_cmd)

    os.system(ray_cmd)

    #print("------- Running damping, mode %d -------" % damp_mode)

    #print(damp_cmd)
    #os.system(damp_cmd)

    # Move back to the working directory
    os.chdir(cwd)

    print('raytracer done')

    return 

# -------------------------------- PARALLELIZE --------------------------------
def parallel_run_rays(time_list, position_list, direction_list, freq_list, output_directory, mds):

    # parallel
    nmbrcores = cpu_count()
    print(nmbrcores, ' run')

    lstarg = zip(time_list, position_list, direction_list, freq_list, output_directory, mds)
    
    with Pool(nmbrcores) as p:
        results = p.starmap(single_run_rays, lstarg)

# -------------------------------- END --------------------------------------