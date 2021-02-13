import numpy as np
import pyglow
from constants_settings import *

"""
pyglow uses: 
:param dn: datetime.datetime object
:param lat: Latitude [degrees]
:param lon: Longitude [degrees]
:param alt: Altitude [km]
"""

class bfieldinfo:
    def __init__(self, time=None, pos=None, B=None, Bvec=None, Bunit=None):
        self.time = None
        self.pos = None
        self.B = None
        self.Bvec = None
        self.Bunit = None
    
    def getBfield(self):
        # we need GEO coordinates sph in km deg deg
        if np.shape(self.pos)[0] > 1:
            Bvec = []
            Bunit = []
            B = []
            for t, pos in zip(self.time, self.pos):
                pt = pyglow.Point(t, pos[1], pos[2], pos[0]-R_E*1e3)
                pt.run_igrf()
                Bvec.append([pt.Bx, pt.By, pt.Bz])
                B.append(pt.B) 
                Bunit.append(np.array([pt.Bx, pt.By, pt.Bz])/pt.B) 
        else:
            pt = pyglow.Point(self.time, self.pos[1], self.pos[2], self.pos[0]-R_E*1e3)
            pt.run_igrf()
            Bvec = [pt.Bx, pt.By, pt.Bz]
            B = pt.B
            Bunit = np.array(Bvec)/B

        # update
        self.Bvec = Bvec
        self.B = B
        self.Bunit = Bunit