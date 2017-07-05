class Box(object):
    """
    Class to define box objects. The attributes are the box length (length),
    a vector containing particle objects (particle) and the number of 
    particles in a box (numParticles)
    """
    def __init__(self, length):
        self.length = length
        self.particle = []
        self.numParticles = 0
