import numpy as np

class Integrators(object):
    """
    """
    def __init__(self, integratorType, timeStep, box):


        self.integratorType = integratorType
        self.timeStep = timeStep
        self.box = box


class VelocityVerlet(Integrators):
    """
    """

    def __init__(self, timeStep, box):
        Integrators.__init__(self, "velocityVerlet", timeStep, box)

    def updatePositions(self):
        """
            
        Parameters
        ----------
	rij2: float

        Returns
        ----------

        Raises
        ----------
        None


        Notes
        ----------
        None

        """

        #print self.box.coordinates
        #raw_input()
        #print self.box.forces
        #raw_input()
        self.box.coordinates = self.box.coordinates \
                + self.box.velocities * self.timeStep \
                + 0.5 * self.box.forces * self.timeStep \
                * self.timeStep #/ self.box.mass

        #print self.box.coordinates
        #raw_input()
        #print self.box.forces
        #raw_input()
        self.box.coordinates = self.box.coordinates - self.box.length * \
                np.round(self.box.coordinates / self.box.length)

    def updateVelocities(self):

        self.box.velocities = self.box.velocities \
     		+ 0.5 * self.box.forces * self.timeStep #/ self.box.mass  
