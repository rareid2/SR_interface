# import needed packages
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Circle
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data

import datetime as dt

from constants_settings import *
from convert_coords import convert2
from raytracer_utils import read_rayfile,readdump, get_yearmiliday 
from bfield import getBline

#---------------------------------------------------------------------------
def rotateplane(rc, tvec_datetime, crs, carsph, units):
    # rotate to be at prime meridian for plotting purposes
    rot_crs = convert2(rc, tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['Re','deg','deg'])

    rot_crs_new = [[rcs[0], rcs[1], 0] for rcs in rot_crs]
    final_crs = convert2(rot_crs_new, tvec_datetime, 'GEO','sph', ['Re','deg','deg'], crs, carsph, units)

    return final_crs
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
def get_Lshells(lshells, tvec_datetime, crs, carsph, units):
    Lshell_fline = []
    for lsh in lshells:
        # needs to be cartesian, output will be in Re
        # plot l shell as if at equator and prime merid
        newt = convert2([[lsh*R_E,0,0]], tvec_datetime, 'GEO', 'sph', ['m','deg','deg'], 'SM', 'car', ['m','m','m'])
        T = getBline(newt[0], tvec_datetime[0])
        
        # repack
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        LshellT = [tvec_datetime[0] for i in range(len(T_repack))]
        
        T_convert = convert2(T_repack, LshellT, 'SM','car',['Re','Re','Re'], crs, carsph, units)
        
        T_repackx = [tt[0] for tt in T_convert]
        T_repackz = [tt[2] for tt in T_convert]
        
        Lshell_fline.append([T_repackx, T_repackz])
    return Lshell_fline
#---------------------------------------------------------------------------

def get_plasmasphere():
    # need to figure out how to change this coordinate system! TODO!
    path2plasma = '/modeldumps/'
    plasma_model_dump = os.path.join(path2plasma, 'model_dump_mode_1_XZ.dat')
    d_xz = readdump(plasma_model_dump)
    Ne_xz = d_xz['Ns'][0, :, :, :].squeeze().T * 1e-6
    Ne_xz[np.isnan(Ne_xz)] = 0

    # Axis spacing depends on how the modeldump was ran
    psize = 10
    px = np.linspace(-10, 10, 200)
    py = np.linspace(-10, 10, 200)

    # Colorbar limits (log space)
    clims = [-2, 5]

    # Plot background plasma (equatorial slice)
    g = plt.pcolormesh(px, py, np.log(Ne_xz), cmap = 'twilight')
    #fig.colorbar(g, ax=ax, orientation="horizontal", pad = 0.2, label= 'Plasmasphere density')

#---------------------------------------------------------------------------
def plotray2D(ray_datenum, raylist, ray_out_dir, crs, carsph, units, plot_kvec=False, show_plot=True):

    # ONLY suports cartesian plotting!

    # convert to desired coordinate system into vector list rays
    ray_coords = []
    k_coords = []
    for r in raylist:
        w = r['w']

        # comes in as SM car in m 
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        tmp_kcoords = list(zip((w/C) * r['n'].x, (w/C) * r['n'].y, (w/C) * r['n'].z))
        
        # convert to a unit vector first
        unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]

        new_kcoords = convert2(unitk, tvec_datetime, 'SM', 'car', ['Re','Re','Re'], crs, carsph, ['Re','Re','Re'])
        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

        # save it
        ray_coords.append(new_coords)
        k_coords.append(new_kcoords)
        
    # -------------------------------- PLOTTING --------------------------------
    # flatten
    ray_coords = ray_coords[0]
    k_coords = k_coords[0]

    ray1 = ray_coords[0]
    ray1_sph = convert2([ray1], tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['m','deg','deg'])

    # find ray starting long
    long = ray1_sph[0][2]
    if long > 180:
        plane_long = 180-(long-180)
        catch = -1
    else:
        plane_long = long
        catch = 1
    try:
        import cartopy
        import cartopy.crs as ccrs
        # make this into a an image (sneaky)
        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=(catch*plane_long)-90, central_latitude=0.0))
        ax.add_feature(cartopy.feature.OCEAN, zorder=0)
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        ax.coastlines()
        plt.savefig(ray_out_dir+'/ccrs_proj.png')
        plt.close()

        img = plt.imread(ray_out_dir+'/ccrs_proj.png', format='png')
        fig, ax = plt.subplots(figsize=(6,6))
        ax.imshow(img, extent=[-1.62,1.62,-1.3,1.3]) # not the most scientific

    except:
        fig, ax = plt.subplots(figsize=(6,6))
        earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
        iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
        ax.add_artist(earth)
        ax.add_artist(iono)

    # plot the rays and kvec
    rotated_rcoords = []
    rotated_kcoords = []
    for rc, kc in zip(ray_coords, k_coords):
        rotated = rotateplane([rc], tvec_datetime, crs, carsph, units)
        rotated_rcoords.append(rotated[0])
        
        rotatedk = rotateplane([kc], tvec_datetime, crs, carsph, units)
        rotated_kcoords.append(rotatedk[0])

    # repack in x,y,z
    rotated_rcoords_x = [rcs[0] for rcs in rotated_rcoords]
    rotated_rcoords_y = [rcs[1] for rcs in rotated_rcoords]
    rotated_rcoords_z = [rcs[2] for rcs in rotated_rcoords]

    # same for k
    kc_x = [kcs[0] for kcs in rotated_kcoords]
    kc_y = [kcs[1] for kcs in rotated_kcoords]
    kc_z = [kcs[2] for kcs in rotated_kcoords]

    if plot_kvec:
        int_plt = len(rotated_rcoords_x)//10
        ax.scatter(rotated_rcoords_x, rotated_rcoords_z, c = 'Black', s = 1, zorder = 103)
        ax.quiver(rotated_rcoords_x[::int_plt], rotated_rcoords_z[::int_plt], -1*catch*np.array(kc_x[::int_plt]), kc_z[::int_plt], color='Black', zorder=104)
        
        for ii,i in enumerate(range(0,len(rotated_rcoords), int_plt)): # not the prettiest but it will work, just trying to annoate
            ax.text(rotated_rcoords_x[i]-0.1, rotated_rcoords_z[i]-0.1, str(ii))
    else:
        ax.scatter(rotated_rcoords_x[0], rotated_rcoords_z[0], c = 'Blue', s=10)
        ax.scatter(rotated_rcoords_x, rotated_rcoords_z, c = 'Black', s = 1, zorder = 103)

    # plot field lines (from IGRF13 model)
    L_shells = [2, 3, 4]  # Field lines to draw
    Lshell_flines = get_Lshells(L_shells, tvec_datetime, crs, carsph, units)
    
    for lfline in Lshell_flines:
        ax.plot(lfline[0], lfline[1], color='Black', linewidth=1, linestyle='dashed')

    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    plt.xlim([0, max(L_shells)])
    plt.ylim([-2, 2])
    plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz')
    plt.savefig(ray_out_dir + '/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_2Dview' + '.png')

    if show_plot:
        plt.show()

    return
# --------------------------------------------------------------------------------------

# ---------------------------------------- STIX PARAM --------------------------------------------
def stix_parameters(ray, t, w):

    B   =  ray['B0'].iloc[t]
    Bmag = np.linalg.norm(B)

    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:])
    Ns   = np.array(ray['Ns'].iloc[t,:])

    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/EPS0/M

    R = 1.0 - np.sum(Wps2/(w*(w + Wcs)))
    L = 1.0 - np.sum(Wps2/(w*(w - Wcs)))
    P = 1.0 - np.sum(Wps2/(w*w))
    S = (R+L)/2.0
    D = (R-L)/2.0

    return R, L, P, S, D
# ---------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
def plotrefractivesurface(ray_datenum, ray, ray_out_dir):
    # set up phi vec
    phi_vec = np.linspace(0,360,int(1e5))*D2R
    w = ray['w']
    
    # comes in as SM car in m 
    tmp_kcoords = list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z))

    # convert to a unit vector first
    unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in ray['time']]
    kcoords = unitk

    int_plt = len(kcoords)//10

    for ti, t in enumerate(ray['time']):
        if ti % int_plt == 0:
            print('plotting a refractive surface')
            # get stix param
            R, L, P, S, D = stix_parameters(ray, ti, w)
            root = -1 # whistler solution

            eta_vec = np.zeros_like(phi_vec)

            # solution from antenna white paper!
            # resangle = np.arctan(np.sqrt(-P/S)) -- not accurate

            # find phi
            B   =  ray['B0'].iloc[ti]
            Bmag = np.linalg.norm(B)
            bunit = B/Bmag
            print(bunit)

            kunit = np.array([kcoords[ti][0], kcoords[ti][1], kcoords[ti][2]])
            print(kunit)

            alpha = np.arccos(np.dot(kunit, bunit))
            alphaedg = float(alpha)*R2D

            alpha_diff_last = 1

            for phi_ind, phi  in enumerate(phi_vec):

                # Solve the cold plasma dispersion relation
                cos2phi = pow(np.cos(phi),2)
                sin2phi = pow(np.sin(phi),2)

                A = S*sin2phi + P*cos2phi
                B = R*L*sin2phi + P*S*(1.0+cos2phi)

                discriminant = B*B - 4.0*A*R*L*P
                n1sq = (B + np.sqrt(discriminant))/(2.0*A)
                n2sq = (B - np.sqrt(discriminant))/(2.0*A)

                # negative refers to the fact that ^^^ B - sqrt
                n1 = np.sqrt(n1sq)
                n2 = np.sqrt(n2sq)

                eta_vec[phi_ind] = n2

                # save the angle near alpha (easier for plotting)
                alpha_diff = np.abs(alpha - phi)
                if alpha_diff < alpha_diff_last:
                    phi_ind_save = phi_ind
                alpha_diff_last = alpha_diff

            # plot it --------------------------
            fig, ax = plt.subplots(1,1, figsize=(6,4))

            # plot the surface
            ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), 'gray', LineWidth = 1, label = 'e + ions')
            
            quiver_mag = eta_vec[phi_ind_save]
            ax.plot([0, quiver_mag*np.sin(float(alpha))], [0, quiver_mag*np.cos(float(alpha))], 'r--')

            xlim1 = -100
            xlim2 = -xlim1
            ylim1 = xlim1
            ylim2 = -ylim1

            ax.set_xlim([xlim1, xlim2])
            ax.set_ylim([ylim1, ylim2])
            ax.set_xlabel('Transverse Refractive Component')
            plt.legend(loc='upper right')
            ax.text(60, 65, r'$\alpha$' +  ' = ' + str(round(alphaedg,2)))

            plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz \n refractive surface ' + str(ti//int_plt))
            plt.savefig(ray_out_dir+'/refractive_surface'+str(ti//int_plt)+'.png')

    return

# ------------------------------------------- END --------------------------------------------

def plot_damp(damplist,thetas, ray_out_dir):

    fig = plt.figure(figsize =(6, 4))

    # damplist is a list of output from dampfiles
    # each dampfile outputs 2 lists, one is a list of ray times
    # the other is a list of ray damps 

    for dn in damplist: # parses each FILE output
        for di, (dtime, ddamp) in enumerate(zip(dn[0],dn[1])): # parses each RAY in file 
            plt.plot(dtime,ddamp,label=str((-1*thetas[di])-180))
 
    """
    # damping attempt
    ray = raylist[0]

    # set up phi vec
    w = ray['w']

    # comes in as SM car in m 
    tmp_kcoords = list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z))

    # convert to a unit vector first
    unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in ray['time']]
    kcoords = unitk

    for ti, t in enumerate(ray['time']):
        # get stix param
        R, L, P, S, D = stix_parameters(ray, ti, w)
        root = -1 # whistler solution

        # find phi
        B   =  ray['B0'].iloc[ti]
        Bmag = np.linalg.norm(B)
        bunit = B/Bmag


        Q    = np.abs(np.array(ray['qs'].iloc[ti,:]))
        M    = np.array(ray['ms'].iloc[ti,:])
        Ns   = np.array(ray['Ns'].iloc[ti,:])

        Wcs   = Q*Bmag/M
        Wps2  = Ns*pow(Q,2)/EPS0/M

        kunit = np.array([kcoords[ti][0], kcoords[ti][1], kcoords[ti][2]])

        alpha = np.arccos(np.dot(kunit, bunit))
        
        alphaedg = float(alpha)*R2D

        # Solve the cold plasma dispersion relation
        cos2phi = pow(np.cos(alpha),2)
        sin2phi = pow(np.sin(alpha),2)

        A = S*sin2phi + P*cos2phi
        B = R*L*sin2phi + P*S*(1.0+cos2phi)

        discriminant = B*B - 4.0*A*R*L*P
        n1sq = (B + np.sqrt(discriminant))/(2.0*A)
        n2sq = (B - np.sqrt(discriminant))/(2.0*A)

        # negative refers to the fact that ^^^ B - sqrt
        n1 = np.sqrt(n1sq)
        n2 = np.sqrt(n2sq)

        eta = n2
        # first coeff
        #An = 1
        #c1 = An / (4*eta*(2*A*eta**2 - B))
        #gamma_constant = 2*np.pi**2 * Wps2[0] / (ray['w'] * kparallel)


        # plot the surface
        #ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), 'gray', LineWidth = 1, label = 'e + ions')
        #quiver_mag = eta_vec[phi_ind_save]
        #ax.plot([0, quiver_mag*np.sin(float(alpha))], [0, quiver_mag*np.cos(float(alpha))], 'r--')

        # i have reasonable confidence in eta, A,B,D,L,P,S,R
    """
    
    plt.xlabel('seconds')
    plt.legend()
    plt.title('geometric factor off')
    plt.savefig(ray_out_dir+'/attempt1_geomoff.png')
    plt.close()

def plotgeomfactor(ray_datenum, raylist, ray_out_dir, crs, carsph, units, show_plot=True):

    # convert to desired coordinate system into vector list rays
    ray_coords = []
    for r in raylist:
        w = r['w']

        # comes in as SM car in m 
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        
        # convert to a unit vector first
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        seconds_passed = [s for s in r['time']]

        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

        # save it
        ray_coords.append(new_coords)
        
    # -------------------------------- PLOTTING --------------------------------
    # flatten
    ray_coords = ray_coords[0]
    r_init = ray_coords[0][0]
    lat_init = ray_coords[0][1]
    
    fig, ax = plt.subplots(figsize=(6,6))
    
    for ti, rc in enumerate(ray_coords):
        # find geometric factor
        gf = r_init * np.cos(lat_init*D2R) / (rc[0] * np.cos(rc[1]*D2R))
        ax.scatter(seconds_passed[ti],gf,color='b')

    plt.xlabel('seconds')
    #plt.gcf().autofmt_xdate()
    plt.ylabel('Geom. Factor')
    plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz')
    plt.savefig(ray_out_dir + '/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_geomfac_RE' + '.png')

    if show_plot:
        plt.show()

    return
# --------------------------------------------------------------------------------------