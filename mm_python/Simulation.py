import numpy as np


class Simulation(object):
    """
    Class containing the main driver for running an MC and MD simulation
    of the Lennard Jones Fluid.
    """
    def __init__(self, method, temperature, steps, printProp, printXYZ,
            ffManager, boxManager, maxDisp = 0.0,
            integrator = None, scaleFreq = 0):
        """
        Constructor of a simulation object.

        Parameters
        ----------
	method: string
	Supported value is "monteCarlo" or "molecularDynamics"

	temperature: float
	System temperature in K

	steps: integer
	Number of Monte Carlo steps or number of Molecular Dynamics
        time steps

	printProp: integer
	Frequency of printing the properties of the system
	(e.g. energy or pressure)

	printXYZ: integer
	Frequency of printing the Cartesian coordinates of the system

	ffManager: ForceFieldManager
        Force Field Manager instance associated to simulation force field

        boxManager: BoxManager
        Box Manager instance associated to simulation box

	maxDisp: float
	Initial maximum displacement of the LJ spheres. Only relevant
        for Monte Carlo simulations.

        integrator: Integrator
        Integrator instance for performing Molecular Dynamics simulations.

        scaleFreq: integer
        Frequency at which atomic velocities will be rescaled to get consistency
        with the target temperature. Relevant only for MD simulations.

        Returns
        ----------
        None

        Raises
        ----------
        None


        Notes
        ----------
        Output energy will be printed in reduced units.
        Maximum displacement will be adjusted to achieve 40% acceptance rate
        in MC simulations.

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
        self.integrator = integrator
        self.scaleFreq = scaleFreq

    def run(self):
        """
        Main driver to perform a Metropolis Monte Carlo simulation of the
        Lennard Jones fuid.
        """
        if self.method == "monteCarlo":
            trajectory = open('trajectory.xyz', 'w')
            self.boxManager.printXYZ(trajectory)
            box = self.boxManager.box
            totalPairEnergy, totalPairVirial = \
                    self.ffManager.getTotalPairEnergyAndVirial(box)
            tailCorrection = self.ffManager.ForceField.getTailCorrection(box)
            pressureCorrection = \
                    self.ffManager.ForceField.getPressureCorrection(box)
            nAccept = 0
            for iStep in range(0, self.steps):
                iParticle = np.random.randint(box.numParticles)
                randomDisplacement = (2.0 * np.random.rand(3) - 1.0)  \
                        * self.maxDisp

                oldPosition = box.coordinates[iParticle].copy()
                oldEnergy, oldVirial = \
                        self.ffManager.getMolPairEnergyAndVirial(iParticle,box)

                box.coordinates[iParticle] += randomDisplacement
                box.coordinates[iParticle] = \
                    box.coordinates[iParticle] - box.length * \
                    np.round(box.coordinates[iParticle]/box.length)

                newEnergy, newVirial = \
                        self.ffManager.getMolPairEnergyAndVirial(iParticle, box)

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
                    totalPairEnergy = totalPairEnergy + dE
                    totalPairVirial = totalPairVirial + (newVirial - oldVirial)
                else:
                    box.coordinates[iParticle] = oldPosition.copy()

                if np.mod(iStep + 1, self.printProp) == 0:
                    accRate = float(nAccept)/(float(iStep)+1) * 100
                    pressure = totalPairVirial
                    pressure += 3.0*box.numParticles/self.beta
                    pressure /= 3.0*np.power(box.length,3)
                    pressure += pressureCorrection

                    totalEnergy = (totalPairEnergy + tailCorrection) \
                            / box.numParticles
                    print(iStep+1, totalEnergy, pressure, accRate, self.maxDisp)

                    if accRate < 38.0:
                        self.maxDisp = self.maxDisp*0.8
                    elif accRate > 42.0:
                        self.maxDisp = self.maxDisp*1.2

                if np.mod(iStep + 1, self.printXYZ) == 0:
                    self.boxManager.printXYZ(trajectory)

            trajectory.close()

        if self.method == "molecularDynamics":

            box = self.boxManager.box

            trajectory = open('trajectory.xyz', 'w')
            self.boxManager.printXYZ(trajectory)

            tailCorrection = self.ffManager.ForceField.getTailCorrection(box)
            pressureCorrection = \
                    self.ffManager.ForceField.getPressureCorrection(box)

            totalPairEnergy, totalPairVirial = \
                    self.ffManager.getTotalPairEnergyAndVirial(box, \
                    populateForces = True)

	    totalEnergy = \
	            (totalPairEnergy + tailCorrection)/ \
	            box.numParticles

            for iStep in range(0,self.steps):

                self.integrator.updatePositions()
                self.integrator.updateVelocities()
                pairEnergy, pairVirial = \
                        self.ffManager.getTotalPairEnergyAndVirial(box, \
                        populateForces = True)
                self.integrator.updateVelocities()

                if np.mod(iStep + 1, self.scaleFreq) == 0:
                    self.boxManager.scaleVelocities(self.temperature)

                if np.mod(iStep + 1, self.printProp) == 0:
                    totalEnergy = (pairEnergy + tailCorrection) \
                            / box.numParticles
                    print(totalEnergy)

                if np.mod(iStep + 1, self.printXYZ) == 0:
                    self.boxManager.printXYZ(trajectory)

            trajectory.close()
