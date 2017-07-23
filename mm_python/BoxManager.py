import numpy as np
import math

class BoxManager(object):
    """
    Class to perform actions on a given box. The main attribute of this
    class is self.box, which represents the box on which actions will
    be performed.
    """
    def __init__(self, box):
        self.box = box

    def addParticles(self, n, method, mass = 0.0):
        """
        Adds n particles to the simulation box using two methods: random
        and lattice. The first one inserts particles randomly. The second
        inserts them in a lattice.

	Parameters
    	----------
    	method : string
	Method to insert particles. Values can be "random" or "lattice".

        n : integer
        Number of molecules to be inserted in the box.

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

        self.box.mass = mass / 6.023e23 * 10.0**-3
        self.box.numParticles = n

        if method == "random":
            self.box.coordinates = (0.5 - np.random.rand(n,3)) * self.box.length

        elif method == "lattice":
            nSide = 1
            self.box.coordinates = np.zeros((n,3))
            while np.power(nSide, 3) < n:
                nSide += 1
            counterPosition = np.zeros((3, ))
            for iParticle in range(0, n):
                self.box.coordinates[iParticle] = \
                    (counterPosition + 0.5)*self.box.length / nSide \
                    - 0.5*self.box.length
                counterPosition[0] += 1
                if counterPosition[0] == nSide:
                    counterPosition[0] = 0
                    counterPosition[1] += 1
                    if counterPosition[1] == nSide:
                        counterPosition[1] = 0
                        counterPosition[2] += 1

            for iParticle in range(0, self.box.numParticles):
                self.box.coordinates[iParticle] -= 0.5


    def assignVelocities(self, temperature):
        """
        This function assigns velocities to atoms according to the
        Maxwell-Boltzmann distribution at a specified temperature.
        It also removes momentum to avoid translational drifts.

	Parameters
    	----------
        temperature: float
        The temperature of the velocity distribution.

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

#        sigma = np.sqrt(1.380648e-23 * temperature / self.box.mass)
#        self.box.velocities = np.random.normal(loc = 0.0, scale = sigma, \
#                size = (self.box.numParticles,3))
#

        self.box.velocities = np.random.rand(self.box.numParticles, 3)

        momentum = np.sum(self.box.mass * self.box.velocities, axis = 0)

        self.box.velocities = self.box.velocities \
            - momentum/(self.box.numParticles * self.box.mass)

    def getConfigFromFile(self, restartFile, mass=0.0):
        """
	Reads in a configuration in xyz format.

	Parameters
    	----------
    	restartFile: string
	Location of the xyz file containing the configuration

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

        self.box.mass = mass / 6.023e23 * 10.0**-3

        with open(restartFile, "r") as f:
            for lineNbr, line in enumerate(f):
                lineSplit = line.split()
                lenLine = len(lineSplit)
                if lenLine == 1:
                    self.box.numParticles = int(lineSplit[0])
                    self.box.coordinates = np.zeros((self.box.numParticles,3))
                if lenLine > 1:
                    iParticle = int(lineSplit[0]) - 1
                    x = float(lineSplit[1])
                    y = float(lineSplit[2])
                    z = float(lineSplit[3])

                    self.box.coordinates[iParticle] = np.array([x,y,z])

    def printXYZ(self, trajectory):
       """
       Prints the current state of the simulation to an xyz file.

       Parameters
       ----------
       toFile: the file where the coordinates will be dumped.

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
       trajectory.write(str(self.box.numParticles)+"\n\n")
       for iParticle in range(0,self.box.numParticles):
           index = "{0:4d}".format(iParticle+1)
           x = "{0:20.15f}".format(self.box.coordinates[iParticle][0])
           y = "{0:20.15f}".format(self.box.coordinates[iParticle][1])
           z = "{0:20.15f}".format(self.box.coordinates[iParticle][2])
           trajectory.write("   ".join([index,x,y,z,"\n"]))

    def scaleVelocities(self, temperature):
       """
       Scales atoms velocities to make them consistent with the
       system temperature.

       Parameters
       ----------
       temperature: float

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
               - self.box.velocities.mean(axis = 0)
       K = 0.5 * np.sum(self.box.velocities * self.box.velocities)
       factor = np.sqrt(1.5 * len(self.box.velocities) * temperature / K)
       self.box.velocities = self.box.velocities * factor
