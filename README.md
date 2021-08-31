# Stanford Raytracer Python Interface

## to set up --> use venv
First, go into SR_interface using:  
```cd SR_interface```  
Next, run:  
```python3 -m venv SR_env``` 

actiavte it (you will need to use this command anytime you're in this folder)  
```source SR_env/bin/activate```    
OR    
```source SR_env/bin/activate.csh```    
(depends on your system)

next to download packages, run the following lines (requirements txt keeps breaking w spacepy):   
```pip3 install numpy```  
```pip3 install matplotlib```  
```pip3 install spacepy```  
```pip3 install requests```  
```pip3 install sgp4```  
```pip3 install pandas```

And a folder for fieldline tracing with IGRF called PyGeopack needs to be installed. To get this, first run:
```export KPDATA_PATH=/path/to/kp```
```export OMNIDATA_PATH=/path/to/omni```
```export GEOPACK_PATH=/path/to/geopack/data```

And now:
```pip3 install PyGeopack --user```  

If the above fails, follow this link and scroll to the 3 steps at the end of the page: 
https://github.com/microsoft/vscode-python/issues/14327

We also need to make sure the libxformd folder is installed correctly. This folder ensures the coordinate transformations are identical to those used in the raytracer. 

First, use the following to update the repo:
```git submodule init```
```git subodule update```

Next, go into the submodule and build it:
```cd libxformd```
```cd xform_double```
```make clean```
```make```
```make shared```
```cd ...``` 

## to run some rays
First, go to 'constants_settings.py' to make sure the settings are correct for your run
You'll need to change the **raytracer_dir** variable to point to the raytracer binary
The run call is:  
```python example_call.py```  

## to run with cartopy   
If you want cartopy (though it will work fine w/o it) visit https://scitools.org.uk/cartopy/docs/latest/ and make sure you have the required dependecies first. 
Next run:  
```pip3 install cartopy``` 

we need to clean up the version of shapely -- it breaks cartopy sometimes  
Run the following lines  
```pip3 uninstall shapely```  
```pip3 install shapely --no-binary shapely```  


## A quick tour around the repo: 
bfield.py -- a quick function to set ray directions, either randomly around the magnetic field line, or at specified angles. If you want field aligned rays, either input (0,0,0) to go north, or use PyGeoPack to go south

constants_settings.py -- a file with all constants and settings (environment and simulation) to be imported into all other files in this repo

convert_coords.py -- a conversion function between coordinate systems that matches the fortran in the ray tracers (uses libxformd), it also supports SM, GEI, and GEO cartesian/spherical

example_call.py -- an example of how to run the raytracer

ray_plots.py -- a huge script with a LOT of plotting tools

raytracer_utils.py -- reads the output of the raytracer and damping

run_model_dump.py -- this sets the plasmasphere (ngo model only) and the dump runs of the model in use
in general you will likely not need this unless you want to plot the plasmasphere

run_rays.py -- call the raytracer! you need Stanford Raytracer built on your computer to use these functions

satellites.py -- class to get satellite positions using most accurate TLEs for the specified date. Also propagates orbits

SLTrack.ini -- settings to access the TLEs for satellites online. Feel free to change to your own account if you want (or if mine is maxing out)
