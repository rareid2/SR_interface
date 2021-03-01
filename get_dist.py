# import needed packages
import numpy as np
import random
from constants_settings import *

# easy function to get a nice distribution defined by # of points and spread
def gaussian_dist(n, spread):

    K = np.array([[1,0], [0,1]])
    m = np.array([0, 0]).reshape(2, 1)
    d = 2
    z = np.random.multivariate_normal(mean=m.reshape(d,), cov=K, size=n)
    y = np.transpose(z)
    dx = y[0] * spread
    dy = y[1] * spread

    return dx, dy


# sample sin distribution for the antenna
def antenna_MC(n):

    thetalist = []
    philist = []

    while len(thetalist) < n:
        xi = random.random()
        th = np.arccos(1-2*xi)
        pxi = random.random()
        if R2D*th > 90:
            continue
        else:
            thetalist.append(R2D*th)
            thetalist.append(R2D*-1*th)

    # generate random azimuth angles uniformly btwn 0 and pi (forward hemi)
    for i in range(int(n)):
        ph = np.pi * random.random()
        philist.append(R2D*ph)

    return thetalist, philist
"""
nrays = 100000
thetas,phis = antenna_MC(nrays)
import matplotlib.pyplot as plt
plt.hist(thetas,bins=50)
plt.ylabel('counts')

plt.xlabel(r'd$\theta$')
plt.show()
"""