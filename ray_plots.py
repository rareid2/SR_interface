# import needed packages
import numpy as np
import datetime as dt
import os 
import seaborn as sns

# plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from inspect import getmembers, isclass

# import funcs from this repo
from constants_settings import *
from convert_coords import convert2
from raytracer_utils import readdump
from bfield import getBline
from satellites import sat


# ------------------------- thank u stack overflow -------------------------
# this function is used when you use plotray2d and turn plotwavenormals on

def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale 
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc
# -----------------------------------------


#------------------------------- rotation for plots --------------------------------------------
def rotateplane(rc, tvec_datetime, crs, carsph, units):
    # rotate to be at prime meridian for plotting purposes
    rot_crs = convert2(rc, tvec_datetime, crs, carsph, units, crs, 'sph', ['Re','deg','deg'])

    rot_crs_new = [[rcs[0], rcs[1], 0] for rcs in rot_crs]
    final_crs = convert2(rot_crs_new, tvec_datetime, crs,'sph', ['Re','deg','deg'], crs, carsph, units)

    return final_crs
#---------------------------------------------------------------------------


# ---------------------------------------- STIX PARAM --------------------------------------------
def stix_parameters(ray, t, w):
    # from austins code

    # get b field
    B   =  ray['B0'].iloc[t]
    Bmag = np.linalg.norm(B)
    # get plasma species
    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:])
    Ns   = np.array(ray['Ns'].iloc[t,:])
    # cyclotron and plasma freq
    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/EPS0/M
    # stix param calculations
    R = 1.0 - np.sum(Wps2/(w*(w + Wcs)))
    L = 1.0 - np.sum(Wps2/(w*(w - Wcs)))
    P = 1.0 - np.sum(Wps2/(w*w))
    S = (R+L)/2.0
    D = (R-L)/2.0

    return R, L, P, S, D
# ---------------------------------------------------------------------------------------------


#-----------------------------  L SHELL TRACING  ----------------------------------------------
def get_Lshells(lshells, tvec_datetime, crs, carsph, units):
    Lshell_fline = []
    for lsh in lshells:
        # needs to be cartesian, output will be in Re
        # plot l shell as if at equator and prime merid
        newt = convert2([[lsh*R_E,0,0]], tvec_datetime, crs, 'sph', ['m','deg','deg'], 'SM', 'car', ['m','m','m'])
        T = getBline(newt[0], tvec_datetime[0],100)
        
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


#---------------------------------------------------------------------------
def plotray2D(ray_datenum, raylist, ray_out_dir, crs, units, md, show_plot=True, return_nmax=False, checklat=None):

    # give it date, the rays, the output directory to save too, the desired coordinates to plot in
    # and the mode

    # show plot will pull up the plot or not

    # return_nmax = true or false, returns the max index of refraction (useful for other stuff)
    # this will make a plot with color indicating the wna
    # if return_nmax is true, you need checklat to be either greater than 0 for northern direction or negative for southern direction

    # !!!!!!! ONLY suports cartesian plotting!

    n_rays_plot = 500 # (dont plot more than 500 on the same plot it looks like chaos)

    # make figures directory if doesnt exist
    if not os.path.exists(ray_out_dir+'/figures'):
        os.makedirs(ray_out_dir+'/figures')

    if return_nmax == True:
        # first, find the max index of refraction (useful for ray spot plots)
        findmax_deg = []
        # convert to desired coordinate system into vector list rays
        ray_coords = []

        # loop through the rays -- NEED ALL THE RAYS
        for ri, r in enumerate(raylist):
            B   =  r['B0'].iloc[0]
            Bmag = np.linalg.norm(B)
            bunit = B/Bmag

            # get the initial wna (already confirmed this is equiv to the initial)
            w = r['w']
            kveci = [(w/C)*r['n'].x.iloc[0],(w/C)*r['n'].y.iloc[0],(w/C)*r['n'].z.iloc[0]]
            kunit = np.array(kveci)/np.linalg.norm(kveci)

            # alpha is wna
            alpha = np.arccos(np.dot(kunit, bunit))
            alphaedg = float(alpha)*R2D
            if checklat < 0:
                alphaedg = 180 - alphaedg
            if alphaedg > 90:
                alphaedg = 90 - (alphaedg - 90)
            findmax_deg.append(alphaedg)

        
            # only want to plot 500 rays to visualize
            if ri < n_rays_plot:
                # comes in as SM car in m 
                tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
                nmag = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
                tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
                new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

                # save it
                ray_coords.append(new_coords) 

        # find max 
        max_wna = max(findmax_deg)  
        max_wna_ind = findmax_deg.index(max(findmax_deg))
        rm = raylist[max_wna_ind]
        nmax = np.sqrt(rm['n'].x.iloc[0]**2 +rm['n'].y.iloc[0]**2+rm['n'].z.iloc[0]**2)

    else:
        # convert to desired coordinate system into vector list rays
        ray_coords = []

        # loop through the rays -- NEED ALL THE RAYS
        for ri, r in enumerate(raylist):
            # only want to plot n_rays_plot rays to visualize -- change this if you want to plot more (but it will get crowded)
            if ri < n_rays_plot:
                # comes in as SM car in m 
                tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
                tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
                new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

                # save it
                ray_coords.append(new_coords) 

    # -------------------------------- PLOTTING --------------------------------
    # just some beuatification
    sns.set_context("notebook", rc={"font.size":16,
                                    "axes.titlesize":18,
                                    "axes.labelsize":12})
    sns.set(font='Franklin Gothic Book',
        rc={'patch.edgecolor': 'w','patch.force_edgecolor': True,})

    # this is all just to get a picture of the Earth on it lol
    ray1 = ray_coords[0][0]
    # this is GEO for the Earth plot
    ray1_sph = convert2([ray1], tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['m','deg','deg'])
    # find ray starting long just for Earth
    long_r = ray1_sph[0][2]
    # convert for cartopy issues
    if long_r > 180:
        plane_long = 180-(long_r-180)
        catch = -1
    else:
        plane_long = long_r
        catch = 1
    try:
        import cartopy
        import cartopy.crs as ccrs
        # make this into a an image (sneaky)

        # chunk of code for the Earth
        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=(catch*plane_long)-90, central_latitude=0.0))
        ax.add_feature(cartopy.feature.OCEAN, zorder=0)
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='w')
        ax.coastlines()
        plt.savefig(ray_out_dir+'/figures/ccrs_proj.png')
        plt.cla()

        import matplotlib.patches as patches
        import matplotlib.cbook as cbook

        img = plt.imread(ray_out_dir+'/figures/ccrs_proj.png', format='png')
        fig, ax = plt.subplots(figsize=(12,12))
        im = ax.imshow(img, extent=[-1.62,1.62,-1.3,1.3],zorder=2) # not the most scientific
        patch = patches.Circle((0, 0), radius=1.02, transform=ax.transData,edgecolor='w')

        im.set_clip_path(patch)

    except: # or just a patch for the Earth
        fig, ax = plt.subplots(figsize=(12,12))
        earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=2)
        iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=1)
        ax.add_artist(earth)
        ax.add_artist(iono)

    # plot the rays
    rotated_rcoords = []
    all_ray_coords = []
    # rotate to be in plane (XZ plane)
    for rayc in ray_coords:
        ray_individ = []
        for rc in rayc:
            rotated = rotateplane([rc], tvec_datetime, crs, carsph, units)
            rotated_rcoords.append(rotated[0])
            ray_individ.append(rotated[0])
        all_ray_coords.append(ray_individ)

    # repack in xz plane
    rotated_rcoords_x = [rcs[0] for rcs in rotated_rcoords]
    rotated_rcoords_z = [rcs[2] for rcs in rotated_rcoords]

    # rays x coords 
    rayx = []
    rayz = []
    for rc in all_ray_coords:
        rx = [float(rr[0]) for rr in rc]
        rz = [float(rr[2]) for rr in rc]
        rayx.append(np.array(rx))
        rayz.append(np.array(rz))

    # plot!
    if return_nmax == True:
        lc = multiline(rayx, rayz, findmax_deg[:n_rays_plot], cmap='coolwarm', lw=2)
        axcb = fig.colorbar(lc)
        axcb.set_label('WNA [deg]')
    else:
        ax.scatter(rotated_rcoords_x, rotated_rcoords_z, c='cornflowerblue', s = 5, zorder = 4)

    # final clean up
    # plot field lines (from IGRF13 model)
    L_shells = [2, 3, 4]  # Field lines to draw
    Lshell_flines = get_Lshells(L_shells, tvec_datetime, crs, carsph, units)
    
    for lfline in Lshell_flines:
        ax.plot(lfline[0], lfline[1], color='dimgrey', linewidth=1, linestyle='dashed',zorder=3)

    # if you wanted to add satellites, here's an example
    """
    # find DSX
    dsx = sat()             # define a satellite object
    dsx.catnmbr = 44344     # provide NORAD ID
    dsx.time = ray_datenum  # set time
    dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag
    vpm = sat()             # define a satellite object    
    vpm.catnmbr = 45120 
    vpm.time = ray_datenum  # set time
    vpm.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

    dsx.propagatefromTLE(sec=0, orbit_dir='future', crs=crs, carsph='car', units=['Re','Re','Re'])
    vpm.propagatefromTLE(sec=0, orbit_dir='future', crs=crs, carsph='car', units=['Re','Re','Re'])

    dsx_spot = dsx.pos
    vpm_spot = vpm.pos
    rotated_dsx = rotateplane(dsx_spot, tvec_datetime, crs, carsph, units)
    rotated_vpm = rotateplane(vpm_spot, tvec_datetime, crs, carsph, units)
    
    plt.scatter(rotated_dsx[0][0],rotated_dsx[0][2],marker='*',c='goldenrod',s=200,zorder=5)
    plt.scatter(rotated_vpm[0][0],rotated_vpm[0][2],marker='*',c='goldenrod',s=200,zorder=6)
    """

    # some clean up 
    plt.xlabel('L (R$_E$)', color='dimgrey')
    plt.ylabel('L (R$_E$)', color='dimgrey')
    #plt.xlim([0, max(L_shells)])
    plt.xlim([0, 4])
    plt.ylim([-1.75, 1.75])

    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['top'].set_color('dimgrey') 
    ax.spines['right'].set_color('dimgrey')
    ax.spines['left'].set_color('dimgrey')
    ax.tick_params(axis='x', colors='dimgrey')    #setting up X-axis tick color to red
    ax.tick_params(axis='y', colors='dimgrey')  #setting up Y-axis tick color to black
    ax.set_facecolor('white')
    
    plt.savefig(ray_out_dir + '/figures/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_2Dview' + str(md) + '.png',bbox_inches='tight')
    rasterize_and_save(ray_out_dir + '/figures/' + dt.datetime.strftime(ray_datenum, '%Y_%m_%d_%H%M%S') + '_' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz' +'_2Dview' + str(md)+'.svg', [ax], fig=fig, dpi=800)

    if show_plot:
        plt.show()
    plt.close()


    if return_nmax:
        return nmax
    else:
        return
# --------------------------------------------------------------------------------------


# ---------------------------------- PLOT REFRACTIVE SURFACE ----------------------------------------------------
def plotrefractivesurface(ray_datenum, raylist, ray_out_dir, ti):

    # ti is the index that you want to plot the refractive surface at
    # ti = 0 is at the ray start point
    # ti = 10 is the tenth time step etc.

    fig, ax = plt.subplots(1,1, figsize=(8,5))
    colors = ['cornflowerblue','g','b']
    for ii,ray in enumerate(raylist):
        # set up phi vec
        phi_vec = np.linspace(0,360,int(1e5))*D2R
        w = ray['w']
        
        # comes in as SM car in m 
        tmp_kcoords = list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z))

        # convert to a unit vector first
        unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
        kcoords = unitk

        # ---------- plot it --------------------------
        print('plotting a refractive surface')

        # get stix param
        R, L, P, S, D = stix_parameters(ray, ti, w)
        root = -1 # whistler solution

        eta_vec = np.zeros_like(phi_vec)

        # solution from antenna white paper!
        resangle = np.arctan(np.sqrt(-P/S)) 

        # find phi
        B   =  ray['B0'].iloc[ti]
        Bmag = np.linalg.norm(B)

        wce = Q_EL*Bmag/M_EL

        #gendrin_angle = 180 - (R2D*np.arccos(2*w/wce))
        #print(gendrin_angle)

        bunit = B/Bmag

        kunit = np.array([kcoords[ti][0], kcoords[ti][1], kcoords[ti][2]])

        # find wna
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


    # plot the surface
    ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), c=colors[ii], LineWidth = 2, label = 'e + ions',zorder=6)
    
    # plot the k-vec
    quiver_mag = eta_vec[phi_ind_save]
    ax.plot([0, quiver_mag*np.sin(float(alpha))], [0, quiver_mag*np.cos(float(alpha))], linestyle=':', color='darkslategrey',zorder=4)
    ax.plot([0, 100*quiver_mag*np.sin(float(resangle))], [0, 100*quiver_mag*np.cos(float(resangle))], linestyle=':', color='lightsteelblue',zorder=4)
    ax.plot([0, -100*quiver_mag*np.sin(float(resangle))], [0, 100*quiver_mag*np.cos(float(resangle))], linestyle=':', color='lightsteelblue',zorder=4)

    # mark the origin
    ax.scatter(0,0,marker='*',s=200,c='lightcoral',zorder=5)

    xlim1 = -750
    xlim2 = -xlim1
    ylim1 = -50
    ylim2 = -ylim1

    # just beautification
    ax.set_xlim([xlim1, xlim2])
    ax.set_ylim([ylim1, ylim2])
    ax.set_xlabel('Transverse Refractive Component',c='dimgrey')
    ax.set_ylabel('Transverse Refractive Component',c='dimgrey')
    ax.set_facecolor('w')
    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['top'].set_color('dimgrey') 
    ax.spines['right'].set_color('dimgrey')
    ax.spines['left'].set_color('dimgrey')
    ax.tick_params(axis='x', colors='dimgrey')    #setting up X-axis tick color to red
    ax.tick_params(axis='y', colors='dimgrey')  #setting up Y-axis tick color to black
    ax.set_ylabel('Transverse Refractive Component',c='dimgrey')
    
    #plt.legend(loc='upper right')
    #ax.text(60, 65, r'$\alpha$' +  ' = ' + str(round(alphaedg,2)))
    
    #plt.title(dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S') + ' ' + str(round(w/(1e3*np.pi*2), 1)) + 'kHz \n refractive surface ')
    #plt.savefig(ray_out_dir+'/refractive_surface'+str(ti)+'.png')
    rasterize_and_save(ray_out_dir + '/figures/' + 'refractive_surface'+str(ti)+'.svg', [ax], fig=fig, dpi=800)

    plt.clf()
    plt.close()


    return

# ------------------------------------------- END --------------------------------------------


# --------------------------- plot geometric factor along ray path! ----------------------------------------------------
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
# ------------------------------------------- END --------------------------------------------


# --------------------------- plot Ne along ray path! ----------------------------------------------------
def plotNe(raylist):
    
    for ri,r in enumerate(raylist):
        Ne = []
        for ti, t in enumerate(r['time']):
            edens = np.array(r['Ns'].iloc[ti])
            Ne.append(np.log10(edens[0]))
        plt.plot(r['time'],Ne,label=['e - model '+str(ri)])

    plt.xlabel('time [sec]')
    plt.xlim([0,0.5])
    plt.ylabel('log(density)')
    plt.show()
    plt.close()

# --------------------------------------------------------------------------------------


# ------------------------------------------- plot plasmasphere density (2D plot) --------------------------------------------
def plot_plasmasphere_2D(mds,kp):
    import cmocean
    import seaborn as sns
    #sns.set(font='Franklin Gothic Book')
    sns.set_context("notebook", rc={"font.size":16,
                                    "axes.titlesize":18,
                                    "axes.labelsize":12})

    ray_out_dir = '/media/rileyannereid/DEMETER/SR_output/2020-08-17_12_00_00'
    if not os.path.exists(ray_out_dir+'/figures'):
        os.makedirs(ray_out_dir+'/figures')

    import cartopy
    import cartopy.crs as ccrs

    # chunk of code for the Earth
    ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0, central_latitude=0.0))
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    ax.coastlines()
    plt.savefig(ray_out_dir+'/figures/ccrs_proj.png')
    plt.cla()

    import matplotlib.patches as patches
    import matplotlib.cbook as cbook

    img = plt.imread(ray_out_dir+'/figures/ccrs_proj.png', format='png')
    
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)#, figsize=(8,12))
    # make this into a an image (sneaky)
    for md,ax in zip(mds,axs):

        im = ax.imshow(img, extent=[-1.75,1.75,-1.62,1.62],zorder=2) # not the most scientific
        patch = patches.Circle((0, 0), radius=1, transform=ax.transData)
        im.set_clip_path(patch)
        
        model_path = 'modeldumps/'
        plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XY_'+str(kp)+'.dat')
        d_xy = readdump(plasma_model_dump)
        plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XZ_'+str(kp)+'.dat')
        d_xz = readdump(plasma_model_dump)
        plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_YZ_'+str(kp)+'.dat')
        d_yz = readdump(plasma_model_dump)
        Ne_xy = d_xy['Ns'][0,:,:,:].squeeze().T
        Ne_xy[np.isnan(Ne_xy)] = 0
        Ne_xz = d_xz['Ns'][0,:,:,:].squeeze().T
        Ne_xz[np.isnan(Ne_xz)] = 0
        Ne_yz = d_yz['Ns'][0,:,:,:].squeeze().T
        Ne_yz[np.isnan(Ne_yz)] = 0

        # Axis spacing depends on how the modeldump was ran
        psize = 10
        px = np.linspace(-10, 10, 1000)
        py = np.linspace(-10, 10, 1000)

        # Plot background plasma 
        g = ax.pcolormesh(px, py, np.log10(Ne_xz*10**(-6)), vmin=0,vmax=4, cmap = cmocean.cm.thermal)

        #earth = plt.Circle((0,0),1,color='k',alpha=1, zorder=100)
        #ax.add_patch(earth)   
        
        ax.set_ylim([-3,3])
        ax.set_xlim([-7,7])
        #ax.set_aspect('equal')
        print('plotted mode ', md)
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.85, 0.11, 0.02, 0.77])
    cb = fig.colorbar(g, cax=cbar_ax,label= 'Plasmasphere Density #/cm^3')
    axs[0].axes.xaxis.set_visible(False)
    axs[0].annotate('Diffusive Eq.',(2,-2.7), c='w')
    axs[1].annotate('GCPM',(4.2,-2.7), c='w')
    axs[1].set_xlabel('L (R$_E$)')

    axs[0].spines['bottom'].set_color('white')
    axs[0].spines['top'].set_color('white') 
    axs[0].spines['right'].set_color('white')
    axs[0].spines['left'].set_color('white')

    axs[1].spines['bottom'].set_color('white')
    axs[1].spines['top'].set_color('white') 
    axs[1].spines['right'].set_color('white')
    axs[1].spines['left'].set_color('white')

    #cbar_ax.spines['bottom'].set_color('white')
    #cbar_ax.spines['top'].set_color('white') 
    #cbar_ax.spines['right'].set_color('white')
    #cbar_ax.spines['left'].set_color('white')
    cb.outline.set_edgecolor('white')


    plt.subplots_adjust(hspace=0.05)
    rasterize_and_save('/home/rileyannereid/workspace/SR_interface/plasmasphere_models.svg', [axs[0],axs[1],cbar_ax], fig=fig, dpi=800)
    plt.close()

# ------------------------------------------- END --------------------------------------------


# ------------------------------------------- plot plasmasphere 1D --------------------------------------------
def plot_plasmasphere_1D(md,kp):
    fig, ax = plt.subplots(1,1,figsize=(4,6))
    # Plot background plasma 
    import seaborn as sns

    # Apply the default theme
    sns.set_theme()

    model_path = 'modeldumps/'
    #for kp in kps:
    #for md in mds:
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XY_'+str(kp)+'.dat')
    d_xy = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_XZ_'+str(kp)+'.dat')
    d_xz = readdump(plasma_model_dump)
    plasma_model_dump = os.path.join(model_path, 'model_dump_mode_'+str(md)+'_YZ_'+str(kp)+'.dat')
    d_yz = readdump(plasma_model_dump)

    labels=  ['e','o+','h']
    cs = ['b','r','g','r']

    for i in range(3):
        print(i)
        Ne_xy = d_xy['Ns'][i,:,:,:].squeeze().T
        Ne_xy[np.isnan(Ne_xy)] = 0
        Ne_xz = d_xz['Ns'][i,:,:,:].squeeze().T
        Ne_xz[np.isnan(Ne_xz)] = 0
        Ne_yz = d_yz['Ns'][i,:,:,:].squeeze().T
        Ne_yz[np.isnan(Ne_yz)] = 0
        # Axis spacing depends on how the modeldump was ran
        psize = 10
        px = np.linspace(-10, 10, 1000)
        py = np.linspace(-10, 10, 1000)

        plt.plot(-px,np.log10(Ne_xz[99,:]),cs[i])

    plt.xlabel('L shell')
    plt.ylabel('Log Equatorial Density (#/cm^3)')

    #plt.xlim([0.5,6])
    #plt.ylim([-1,12])
    #plt.legend(['ngo','diff eq','gcpm'])
    #plt.legend([0,4])
    plt.legend(labels)
    plt.grid()
    plt.show()
    plt.close()

# ------------------------------------------- END --------------------------------------------


def plot_density_alongpath(ray_datenum, ray, ray_out_dir,t_save):
    fig = plt.figure(figsize=(3.5,5))
    xs = np.linspace(0,1,len(ray['time']))
    #print(ray['time'],len(ray['time']))
    for ti,t in enumerate(ray['time']):
        raydens = ray['Ns'].iloc[ti]
        Ne = raydens[28]
        Nps = raydens[29]
        Nhe = raydens[30]
        Nos = raydens[31]

        plt.scatter(xs[ti],np.log10(Ne),c='k',s=10)
        plt.scatter(xs[ti],np.log10(Nps),c='y',s=10)
        plt.scatter(xs[ti],np.log10(Nhe),c='g',s=10)
        plt.scatter(xs[ti],np.log10(Nos),c='b',s=10)

        if ti==t_save:
            plt.scatter(xs[ti],np.log10(Ne),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nps),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nhe),c='r',zorder=100)
            plt.scatter(xs[ti],np.log10(Nos),c='r',zorder=100)

    plt.ylabel('log density #/m^3')
    #plt.legend(['electrons (protons)', 'helium ions', 'oxygen ions'])
    #plt.ylim([6,9.75])
    #plt.show()
    plt.savefig(ray_out_dir+'/density_path'+str(t_save)+'.png')
    plt.close()

def calc_fce(ray, t):

    # get b field
    B   =  ray['B0'].iloc[t]
    Bmag = np.linalg.norm(B)
    # get plasma species
    # only want the electron charge and mass
    Q    = np.abs(np.array(ray['qs'].iloc[t,0]))
    M    = np.array(ray['ms'].iloc[t,0])

    # cyclotron freq
    Wce   = Q*Bmag/M

    # convert to Hz
    fce = Wce / (2*np.pi)
    return fce

def calc_gendrin(ray,t):

    f = ray['w'] / (2 * np.pi)
    fce = calc_fce(ray,t)

    th_g = 180 - np.arccos(f/fce)

    th_g = th_g * R2D

    return th_g

def rasterize_and_save(fname, rasterize_list=None, fig=None, dpi=None,
                       savefig_kw={}):
    """Save a figure with raster and vector components
    This function lets you specify which objects to rasterize at the export
    stage, rather than within each plotting call. Rasterizing certain
    components of a complex figure can significantly reduce file size.
    Inputs
    ------
    fname : str
        Output filename with extension
    rasterize_list : list (or object)
        List of objects to rasterize (or a single object to rasterize)
    fig : matplotlib figure object
        Defaults to current figure
    dpi : int
        Resolution (dots per inch) for rasterizing
    savefig_kw : dict
        Extra keywords to pass to matplotlib.pyplot.savefig
    If rasterize_list is not specified, then all contour, pcolor, and
    collects objects (e.g., ``scatter, fill_between`` etc) will be
    rasterized
    Note: does not work correctly with round=True in Basemap
    Example
    -------
    Rasterize the contour, pcolor, and scatter plots, but not the line
    >>> import matplotlib.pyplot as plt
    >>> from numpy.random import random
    >>> X, Y, Z = random((9, 9)), random((9, 9)), random((9, 9))
    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2)
    >>> cax1 = ax1.contourf(Z)
    >>> cax2 = ax2.scatter(X, Y, s=Z)
    >>> cax3 = ax3.pcolormesh(Z)
    >>> cax4 = ax4.plot(Z[:, 0])
    >>> rasterize_list = [cax1, cax2, cax3]
    >>> rasterize_and_save('out.svg', rasterize_list, fig=fig, dpi=300)
    """

    # Behave like pyplot and act on current figure if no figure is specified
    fig = plt.gcf() if fig is None else fig

    # Need to set_rasterization_zorder in order for rasterizing to work
    zorder = -5  # Somewhat arbitrary, just ensuring less than 0

    if rasterize_list is None:
        # Have a guess at stuff that should be rasterised
        types_to_raster = ['QuadMesh', 'Contour', 'collections']
        rasterize_list = []

        print("""
        No rasterize_list specified, so the following objects will
        be rasterized: """)
        # Get all axes, and then get objects within axes
        for ax in fig.get_axes():
            for item in ax.get_children():
                if any(x in str(item) for x in types_to_raster):
                    rasterize_list.append(item)
        print('\n'.join([str(x) for x in rasterize_list]))
    else:
        # Allow rasterize_list to be input as an object to rasterize
        if type(rasterize_list) != list:
            rasterize_list = [rasterize_list]

    for item in rasterize_list:

        # Whether or not plot is a contour plot is important
        is_contour = (isinstance(item, matplotlib.contour.QuadContourSet) or
                      isinstance(item, matplotlib.tri.TriContourSet))

        # Whether or not collection of lines
        # This is commented as we seldom want to rasterize lines
        # is_lines = isinstance(item, matplotlib.collections.LineCollection)

        # Whether or not current item is list of patches
        all_patch_types = tuple(
            x[1] for x in getmembers(matplotlib.patches, isclass))
        try:
            is_patch_list = isinstance(item[0], all_patch_types)
        except TypeError:
            is_patch_list = False

        # Convert to rasterized mode and then change zorder properties
        if is_contour:
            curr_ax = item.ax.axes
            curr_ax.set_rasterization_zorder(zorder)
            # For contour plots, need to set each part of the contour
            # collection individually
            for contour_level in item.collections:
                contour_level.set_zorder(zorder - 1)
                contour_level.set_rasterized(True)
        elif is_patch_list:
            # For list of patches, need to set zorder for each patch
            for patch in item:
                curr_ax = patch.axes
                curr_ax.set_rasterization_zorder(zorder)
                patch.set_zorder(zorder - 1)
                patch.set_rasterized(True)
        else:
            # For all other objects, we can just do it all at once
            curr_ax = item.axes
            curr_ax.set_rasterization_zorder(zorder)
            item.set_rasterized(True)
            item.set_zorder(zorder - 1)

    # dpi is a savefig keyword argument, but treat it as special since it is
    # important to this function
    if dpi is not None:
        savefig_kw['dpi'] = dpi

    # Save resulting figure
    fig.savefig(fname, **savefig_kw, tight_layout=True)