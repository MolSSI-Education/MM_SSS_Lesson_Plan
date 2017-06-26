import numpy as np


class particle(object):
    def __init__(self):

        self.position = np.zeros((3, ), dtype=np.float64)
        self.parms = np.empty((2, ), dtype=np.float64)
        self.energy = 0.0
