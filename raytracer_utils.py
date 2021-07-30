# import packages needed
import numpy as np
import pandas as pd
import os

def get_yearmiliday(ray_datenum):
    # convert for raytracer settings for convenience
    days_in_the_year = ray_datenum.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(ray_datenum.year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = ray_datenum.hour*3.6e6 + ray_datenum.minute*6.0e4 + ray_datenum.second*1.0e3

    return yearday, milliseconds_day


def read_rayfiles(directory, freq, latmin, latmax, lonmin, lonmax):
    '''read rayfiles from a directory '''

    # if (lonmin_in < 0):
    #     lonmin = 360 + lonmin_in
    # else:
    #     lonmin = lonmin_in
    # if (lonmax_in < 0):
    #     lonmax = 360 + lonmax_in
    # else:
    #     lonmax = lonmax_in

    # if (lonmin > lonmax):
    #     tmp = lonmin
    #     lonmin = lonmax
    #     lonmax = tmp

    print("lonmin: ", lonmin)
    print("lonmax: ", lonmax)
    print("hello")
    out = []
    for root, dirs, files in os.walk(os.path.join(directory, "f_%d"%freq)):
        for f in files:
            if f.startswith('ray') and f.endswith('.ray'):
                row = (f.split(".")[0]).split("_")
                print(row)
                cfreq = float(row[1])
                clat  = float(row[2])
                clon  = float(row[3])
                
                if ( (cfreq == freq) and (clat >= latmin) and (clat <= latmax) and
                (clon >= lonmin) and (clon <= lonmax) ):
                    print("reading rayfile at", clat, clon, f)
                    tmp = read_rayfile(os.path.join(root, f))

                    # check damping:
                    dpath = os.path.join(root, 'damp_%d_%d_%d.ray'%(cfreq, clat, clon))

                    if os.path.exists(dpath):
                        dtmp = read_damp(dpath)

                        for ind in range(len(tmp)):
                            tmp[ind]['damping'] = dtmp[ind]['damping']
                    
                    out.extend(tmp)
                    # out.extend(read_rayfile(os.path.join(root,f)))

    return out


def read_damp(dampfile):
    x = pd.read_table(dampfile,delim_whitespace=True,header=None,names=['ray','co','no'],engine='python')

    raynums = np.unique(x[0])
    numrays = len(raynums)

    out = []

    for ii in range(numrays):
        tmp = x[x[0]==raynums[ii]]
        tmp.reset_index(drop=True)
        data = dict()
        data['time'] = tmp[1]
        data['damping'] = tmp[2]
        out.append(data)

    return out

def read_damp_simple(dampfile):
    f = open(dampfile)
    raynums = []
    lines = []
    for line in f:
        lines.append(line)
        sline = line.split()
        raynum = sline[0]
        raynums.append(raynum)
    f.close()


    # find number of rays
    # converting our list to set
    new_set = set(raynums)
    print(len(new_set), ' rays found')
    rays_t = []
    rays_d = []
    for rn in range(1,len(new_set)+1):
        raytime = []
        raydamp = []
        for line in lines:
            sline = line.split()
            if int(sline[0]) == rn:
                raytime.append(float(sline[1]))
                raydamp.append(float(sline[2]))
        rays_t.append(raytime)
        rays_d.append(raydamp)

    return [rays_t, rays_d]

def read_rayfile(rayfile):
    ''' Load output from Forest's raytracer'''
    x = pd.read_csv(rayfile,delim_whitespace=True, header=None)
    # % fields: raynum, stopcond, time, pos(3), vprel(3), vgrel(3), n(3),
    # % B0(3), Nspec, qs(Nspec), ms(Nspec), Ns(Nspec), nus(Nspec)

    raynums = np.unique(x[0])
    numrays = len(raynums)

    out = []
    for ii in range(numrays):
        tmp = x[x[0]==raynums[ii]]
        tmp.reset_index(drop=True)
        #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        data = dict()
        data['time'] = tmp[2]
        data['pos'] = tmp.loc[:,3:5]
        data['pos'].columns=['x','y','z']
        data['vprel'] = tmp.loc[:,6:8]
        data['vprel'].columns=['x','y','z']
        data['vgrel'] = tmp.loc[:,9:11]
        data['vgrel'].columns=['x','y','z']
        data['n'] = tmp.loc[:,12:14]
        data['n'].columns=['x','y','z']
        data['B0'] = tmp.loc[:,15:17]
        data['B0'].columns=['x','y','z']
        
        # #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        # data = dict()
        # data['time'] = tmp[2]

        # data['pos']   = coord.Coords(zip(tmp.loc[:,3], tmp.loc[:,4],  tmp.loc[:,5]), 'SM','car')

        # data['pos']   = data['pos'].
        # data['vprel'] = tmp.loc[:,6:8]
        # data['vprel'].columns=['x','y','z']
        # data['vgrel'] = tmp.loc[:,9:11]
        # data['vgrel'].columns=['x','y','z']
        # data['n'] = tmp.loc[:,12:14]
        # data['n'].columns=['x','y','z']
        # data['B0'] = tmp.loc[:,15:17]
        # data['B0'].columns=['x','y','z']

        data['w'] = tmp.iloc[0,18]       # Frequency
        data['Nspec'] = tmp.iloc[0,19]   # Number of species
        data['stopcond'] = tmp.iloc[0,1] # Stop condition
        data['qs'] =  tmp.loc[:,20 + 0*data['Nspec']:20 + 1*data['Nspec'] - 1]
        data['ms'] =  tmp.loc[:,20 + 1*data['Nspec']:20 + 2*data['Nspec'] - 1]
        data['Ns'] =  tmp.loc[:,20 + 2*data['Nspec']:20 + 3*data['Nspec'] - 1]
        data['nus']=  tmp.loc[:,20 + 3*data['Nspec']:20 + 4*data['Nspec'] - 1]

        # read damping data too, if we have it:
        if (tmp.shape[1] > 20+4*data['Nspec']):
            print('we have it')
            data['damping'] = tmp.loc[:,20+4*data['Nspec']]

            
        out.append(data)
    return out

def readdump(filename):
    ''' Reads the dump files from Forests raytracer
        % x - vector of x coordinates in meters
        % y - vector of y coordinates in meters
        % z - vector of z coordinates in meters
        % 
        % qs - array(ispecies, ix, iy, iz) of charges for each species
        % Ns - array(ispecies, ix, iy, iz) of number densities in m^-3 
        % Ms - array(ispecies, ix, iy, iz) of masses in kg
        % nus - array(ispecies, ix, iy, iz) of collision frequencies in s^-1
        % B0 - array(icomponent, ix,iy, iz) containing the background magnetic
    '''
    out = dict()
    
    f = open(filename,'r')
    line = f.readline().split()
    nspec = int(line[0])
    nx = int(line[1])
    ny = int(line[2])
    nz = int(line[3])
    
    line = f.readline().split()
    minx = float(line[0])
    maxx = float(line[1])
    miny = float(line[2])
    maxy = float(line[3])
    minz = float(line[4])
    maxz = float(line[5])
    
#     print maxz
    
    dat = [float(x) for x in f.readlines()]
#     print np.shape(dat)
#     print dat[0:10]
    dat = np.reshape(dat,[nspec*4 + 3, nx, ny, nz],order='f')
#     print np.shape(dat)
    
    out['qs'] = dat[0*nspec:(1*nspec),:,:,:]
    out['Ns'] = dat[1*nspec:(2*nspec),:,:,:]
    out['Ms'] = dat[2*nspec:(3*nspec),:,:,:]
    out['nus']= dat[3*nspec:(4*nspec),:,:,:]
    out['B0'] = dat[4*nspec:-1,:,:,:]
    
    out['x'] = np.linspace(minx, maxx, nx)
    out['y'] = np.linspace(miny, maxy, ny)
    out['z'] = np.linspace(minz, maxz, nz)
    
    return out

def read_input_jobs(fname):
    f = open(fname)
    lats = []
    lons = []
    psds = []
    for line in f:
        lines = line.split()
        lons.append(float(lines[0]))
        lats.append(float(lines[1]))
        psds.append(float(lines[2]))
    return lats, lons, psds

def read_damp_matlab(fname):
    f = open(fname)
    damp_data = []
    for line in f:
        lines = line.split()
        if lines[0] == 't':
            continue
        else:
            damp_data.append(float(lines[1]))
    return damp_data

# for really large ray files (several gigs)
def read_bigrayfile(rayfile):
    ''' Load output from Forest's raytracer'''
    num=0

    ray_data = []
    iray_data = []
    
    last_raynum = 1
    last_line = 0
    bad_ray = 0
    with open(rayfile) as a_file:
        line_count = 0
        for line in a_file:
            lines = line.split()
            ray_num = int(lines[0])

            if int(lines[1])!=1:
                if ray_num == last_raynum:
                    pass
                else:
                    bad_ray+=1

            if line_count == 0: # get the first ray
                iray_data.append([float(lines[3]),float(lines[4]),float(lines[5])])
            line_count+=1
            
            if ray_num == last_raynum:
                pass
            elif ray_num != last_raynum:
                ray_data.append([float(last_line[3]),float(last_line[4]),float(last_line[5])])
                iray_data.append([float(lines[3]),float(lines[4]),float(lines[5])])

            last_line = lines
            last_raynum = ray_num
                
        # get the last ray!
        ray_data.append([float(lines[3]),float(lines[4]),float(lines[5])])
        #print('imma keep goin')
        #pass

    print(bad_ray, ray_num)
    a_file.close()
    return ray_data, iray_data

def read_bigrayfile_in(rayfile):
    ''' Load output from Forest's raytracer'''
    num=0

    ray_data = []
    last_raynum = 0 # catch the first ray
    last_line = 0
    bad_ray = 0
    with open(rayfile) as a_file:

        for line in a_file:
            lines = line.split()
            ray_num = int(lines[0])

            if int(lines[1])!=1:
                if ray_num == last_raynum:
                    pass
                else:
                    bad_ray+=1

            if ray_num == last_raynum:
                pass
            elif ray_num != last_raynum: # save initial position!
                ray_data.append([float(lines[3]),float(lines[4]),float(lines[5])])

            last_line = lines
            last_raynum = ray_num
                
        # get the last ray!
        #ray_data.append([float(lines[3]),float(lines[4]),float(lines[5])])
        #print('imma keep goin')
        #pass

    print(bad_ray, ray_num)
    a_file.close()
    return ray_data