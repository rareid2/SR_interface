# Stanford Raytracer Python Interface

## to set up --> use venv
run:
```python3 -m venv SR_env``` 

actiavte it (you will need to use this command anytime you're in this folder)
```source SR_env/bin/activate``` 

next to download packages, run: 
```pip install -r requirements.txt``` 

we need to clean up the version of shapely -- it breaks cartopy sometimes
Run the following lines
```pip uninstall shapely```
```pip install shapely --no-binary shapely```

## to run some rays
First, go to 'constants_settings.py' to make sure the settings are correct for your run
You'll need to change the **raytracer_dir** variable to point to the raytracer binary

Next, go to 'example_call.py' 
Here you'll find info on creating a 'ray_inpfile.txt' that the raytracer will run
Output will be stored as specified in the file

You can also plot the ray after it runs using 'ray_plots.py'