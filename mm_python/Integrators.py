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
    This implements the integration routines for the velocity
    verlet method. Two functions are included, update positions
    and update velocities. 
    """

    def __init__(self, timeStep, box):
        Integrators.__init__(self, "velocityVerlet", timeStep, box)

    def updatePositions(self):
        """
        Updates the positions according to the Velocity Verlet method.
        It also applies PBC.

        Parameters
        ----------
        None

        Returns
        ----------
        None

        Raises
        ----------
        None


        Notes
        ----------
        None

        """

        self.box.coordinates = self.box.coordinates \
                + self.box.velocities * self.timeStep \
                + 0.5 * self.box.forces * self.timeStep \
                * self.timeStep 

        self.box.coordinates = self.box.coordinates - self.box.length * \
                np.round(self.box.coordinates / self.box.length)

    def updateVelocities(self):
        """
        Updates the velocities according to the Velocity Verlet method
        using the current forces.

        Parameters
        ----------
        None

        Returns
        ----------
        None

        Raises
        ----------
        None


        Notes
        ----------
        None

        """

        self.box.velocities = self.box.velocities \
                + 0.5 * self.box.forces * self.timeStep 
