# import needed packages
import numpy as np

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

    for i in range(n):
        xi = random.random()
        th = np.arccos(1-2*xi)
        pxi = random.random()
        thetalist.append(R2D*th)

    # generate random azimuth angles uniformly btwn 0 and pi (forward hemi)
    for i in range(int(n)):
        ph = np.pi * random.random()
        philist.append(R2D*ph)

    return thetalist, philist