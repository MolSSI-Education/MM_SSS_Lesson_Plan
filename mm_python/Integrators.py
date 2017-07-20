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

        self.box.coordinates = self.box.coordinates \
                + self.box.velocities * self.timeStep * 10.0**10 \
                + 0.5 * self.box.forces * 1.380648e-23 * self.timeStep \
                * self.timeStep * 10.0**20\
                / self.box.mass

        self.box.coordinates = self.box.coordinates - self.box.length * \
                np.round(self.box.coordinates / self.box.length)

    def updateVelocities(self):

        self.box.velocities = self.box.velocities \
     		+ 0.5 * self.box.forces * 1.380648e-23 * self.timeStep \
                * 10.0**10 / self.box.mass  
