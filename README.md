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
```pip3 install PyGeopack```  

We also need to make sure the libxformd folder is installed correctly. This folder ensures the coordinate transformations are identical to those used in the raytracer. Use the following:
```cd libxformd```
```cd xform_double```
```make clean```
```make```
```make shared```

Next, make sure you have the folder modeldumps  
```mkdir modeldumps```  

## to run some rays
First, go to 'constants_settings.py' to make sure the settings are correct for your run
You'll need to change the **raytracer_dir** variable to point to the raytracer binary
The run call is:  
```example_call.py```  

## to run with cartopy   
If you want cartopy (this will work fine w/o it) visit https://scitools.org.uk/cartopy/docs/latest/ and make sure you have the required dependecies first. 
Next run:  
```pip3 install cartopy``` 

we need to clean up the version of shapely -- it breaks cartopy sometimes  
Run the following lines  
```pip3 uninstall shapely```  
```pip3 install shapely --no-binary shapely```  