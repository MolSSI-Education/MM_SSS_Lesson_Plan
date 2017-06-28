import numpy as np


class Particle(object):
    def __init__(self):

        self.position = np.zeros((3, ), dtype=np.float64)
        self.parms = np.empty((2, ), dtype=np.float64)
