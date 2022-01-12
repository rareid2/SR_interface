# import packages needed to define settings
import numpy as np

# Constants (matching the original Fortran code)
EPS0 = 8.854187817e-12        # C^2/Nm^2 
MU0 = np.pi * 4e-7              
C = np.sqrt(1.0/(EPS0*MU0))   # m/s
R_E = 6371.2e3                # m 
D2R = np.pi / 180.0           # rad
R2D = 180.0 / np.pi           # deg
H_IONO = 1000e3               # m
Hz2Rad = 2.*np.pi             # rad
Rad2Hz = 1./Hz2Rad            # Hz
Q_EL = 1.602e-19              # C
M_EL = 9.1e-31                # kg
M_P = 1.67e-27                # kg

# where is the raytracer?? point to bin directory
raytracer_dir = '/home/rileyannereid/workspace/Stanford_Raytracer/bin'

#  --------------------------- CHANGE ENV SETTINGS HERE  --------------------------
# Environmental parameters
Kp = 1
AE = 1.6
Pdyn = 4
Dst = 1.0
ByIMF = 0.0
BzIMF = -5
# Tsyganenko correction parameters
W = [0.132, 0.303, 0.083, 0.070, 0.211, 0.308]  # Doesn't matter if we're not using Tsyg

# Simulation parameters
t_max = 10        # Maximum duration in seconds
dt0 = 1e-3       # Initial timestep in seconds
dtmax = 0.1      # Maximum allowable timestep in seconds
root = 2         # Which root of the Appleton-Hartree equation
                 # (1 = negative, 2 = positive)
                 # (2 = whistler in magnetosphere)
fixedstep = 0    # Don't use fixed step sizes, that's a bad idea.
maxerr = 1.0e-4  # Error bound for adaptive timestepping
maxsteps = 1e3   # Max number of timesteps (abort if reached)
use_IGRF = 1     # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 1     # Use the Tsyganenko magnetic field model corrections
minalt = R_E+500e3   # cutoff altitude in meters

# Should we include a geometric focusing term in the damping?
include_geom_factor = 0  # 1 for yes

# interpolation parameters for mode 4
scattered_interp_window_scale = 1.2
scattered_interp_order = 2
scattered_interp_exact = 0  # Try 0 if there's weirdness at discontinuities
scattered_interp_local_window_scale = 5

# ngo config
configfile = 'ngoconfig.in'

#  --------------------------- END CHANGE ENV SETTINGS HERE  --------------------------