import numpy as np
from .Particle import Particle

class BoxManager(object):
    """
    Class to perform actions on a given box. The main attribute of this
    class is self.box, which represents the box on which actions will
    be performed. 
    """
    def __init__(self, box):
        self.box = box

    def addParticles(self, n, method):
        """
        Adds n particles to the simulation box using two methods: random
        and lattice. The first one inserts particles randomly. The second
        inserts them in a lattice. 

	Parameters
    	----------
    	method : string
	Method to insert particles. Values can be "random" or "lattice".

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

        self.box.numParticles = n 

        if method == "random":
            #for iParticle in range(0, n):   
                #self.box.particle.append(Particle())
                #newPosition = (0.5 - np.random.rand(3)) * self.box.length
                #self.box.particle[iParticle].position = newPosition
                #self.box.numParticles = self.box.numParticles + 1
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

    def getConfigFromFile(self, restartFile):
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
           x = "{0:10.5f}".format(self.box.coordinates[iParticle][0])
           y = "{0:10.5f}".format(self.box.coordinates[iParticle][1])
           z = "{0:10.5f}".format(self.box.coordinates[iParticle][2])
           trajectory.write("   ".join(["Ar",x,y,z,"\n"]))
