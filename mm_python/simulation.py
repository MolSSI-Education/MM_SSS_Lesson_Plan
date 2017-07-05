import numpy as np

kB = 0.0083144621    #Boltzmann constant, kJ/molK

class Simulation(object):
    """
    Class containing the main driver for running an MC simulation
    of the Lennard Jones Fluid.
    """
    def __init__(self, method, temperature, steps, printProp, printXYZ, maxDisp, ffManager, boxManager):
        """
        Constructor of a simulation object.
 
        Parameters
        ----------
	method: string
	Only supported value is "monteCarlo".

	temperature: float
	System temperature in K

	steps: integer
	Number of Monte Carlo steps

	printProp: integer
	Frequency of printing the thermodynamic properties of the system
	(e.g. energy)

	printXYZ: integer
	Frequency of printing the Cartesian coordinates of the system

	maxDisp: float
	Initial maximum displacement of the LJ spheres

	ffManager: ForceFieldManager
        Force Field Manager instance associated to simulation force field

        boxManager: BoxManager
        Box Manager instance associated to simulation box
	

        Returns
        ----------
        None

        Raises
        ----------
        None


        Notes
        ----------
        Output energy will be printed in reduced units.
        Maximum displacement will be adjusted to achieve 40% acceptance rate.


        """


        self.method = method
        self.temperature = temperature
        self.steps = steps
        self.printProp = printProp
        self.printXYZ = printXYZ
        self.maxDisp = maxDisp
        self.ffManager = ffManager
        self.boxManager = boxManager
        self.beta = 1.0 / (self.temperature)

    def run(self):
        """
        Main driver to perform a Metropolis Monte Carlo simulation of the
        Lennard Jones fuid.
        """
        if self.method == "monteCarlo":
            box = self.boxManager.box
            pairEnergy = self.ffManager.getPairEnergy(box)
            tailCorrection = self.ffManager.ForceField.getTailCorrection(box)
            nAccept = 0
            for iStep in range(0, self.steps):
                iParticle = np.random.randint(box.numParticles)
                randomDisplacement = (2.0 * np.random.rand(3) - 1.0) * \
                        self.maxDisp

                oldPosition = box.particle[iParticle].position.copy()
                oldEnergy = self.ffManager.getMolEnergy(iParticle,box)
                
                box.particle[iParticle].position += randomDisplacement
                box.particle[iParticle].position = \
                    box.particle[iParticle].position - box.length * \
                    np.round(box.particle[iParticle].position/box.length)

                newEnergy = self.ffManager.getMolEnergy(iParticle, box)
                dE = newEnergy - oldEnergy

                accept = False
                if dE <= 0.0:
                    accept = True
                else:
                    randomNumber = np.random.rand(1)[0]
                    factor = self.beta * dE
                    pAcc = np.exp(-factor)
                    if randomNumber < pAcc:
                        accept = True

                if accept:
                    nAccept = nAccept + 1
                    pairEnergy = pairEnergy + dE
                else:
                    box.particle[iParticle].position = oldPosition.copy()

                if np.mod(iStep + 1, self.printProp) == 0:
                    accRate = float(nAccept)/(float(iStep)+1) * 100
                    totalEnergy = \
                            (pairEnergy + tailCorrection)/ \
                            (self.ffManager.ForceField.parms[1]* \
                            box.numParticles)
                    print(iStep+1, totalEnergy, accRate, self.maxDisp)

                    if accRate < 38.0:
                        self.maxDisp = self.maxDisp*0.8
                    elif accRate > 42.0:
                        self.maxDisp = self.maxDisp*1.2

                if np.mod(iStep, self.printXYZ) == 0:
                    pass

