import numpy as np

class Box(object):
    """
    Class to define box objects. The attributes are the box length (length),
    velocities and coordinates.
    """
    def __init__(self, length):
        self.length = length
        self.velocities = None
        self.coordinates = None
