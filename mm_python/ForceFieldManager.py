import numpy as np


class ForceFieldManager(object):
    def __init__(self, ForceField):
        self.ForceField = ForceField


    def getMolEnergy(self, iParticle, box):
        """
	Computes the inter-molecular interaction of a given particle
	in a given box.

        Parameters
        ----------
	box: box
	The box containing the particle.

	iParticle: integer
	Index of the particle.

        Returns
        ----------
        eij: float
	Inter-molecular energy of particle iParticle in box box


        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        iPosition = box.coordinates[iParticle]
        eij = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.coordinates[jParticle]
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.ForceField.cutoff2:
                eij += self.ForceField.evaluate(rij2)
        return eij



    def getMolPairEnergyAndVirial(self, iParticle, box, populateForces = False):

        iPosition = box.coordinates[iParticle]
        ePair = 0.0
        wPair = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.coordinates[jParticle]
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.ForceField.cutoff2:
                ePair += self.ForceField.evaluate(rij2)
                wPair += self.ForceField.getPairVirial(rij2)
                if populateForces == True:
                    box.forces[iParticle] += wPair * rij / rij2

        return ePair, wPair



    def getTotalPairEnergyAndVirial(self, box, populateForces = False):
        """
	Computes the inter-molecular energy of a box.

        Parameters
        ----------
	box: box
	The box containing the particles.

        Returns
        ----------
        ePair: float
	Total inter-molecular energy of box.

        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        if populateForces == True:
            box.forces = np.zeros((box.numParticles,3))

        ePair = 0.0
        wPair = 0.0
        
        for iParticle in range(0, box.numParticles):
            eInter, wInter = self.getMolPairEnergyAndVirial \
                (iParticle, box, populateForces)
            wPair = wPair + wInter
            ePair = ePair + eInter

        wPair = wPair/2.0
        ePair = ePair/2.0
        return ePair, wPair

        #if populateForces == True:
        #    box.forces = np.zeros((box.numParticles,3))

        #ePair = 0.0
        #wPair = 0.0
        #

        #for iParticle in range(0, box.numParticles):
        #    eInter = self.getMolEnergy(iParticle, box)
        #    wInter = self.getMolVirial(iParticle, box, populateForces)
        #    wPair = wPair + wInter
        #    ePair = ePair + eInter

        #wPair = wPair/2.0
        #ePair = ePair/2.0
        #return ePair, wPair

    def getMolVirial(self, iParticle, box, populateForces = False):
        """
	Computes the virial interaction of a given particle

        Parameters
        ----------
	box: box
	The box containing the particle.

	iParticle: integer
	Index of the particle.

        Returns
        ----------
        wPair: float
        Computed virial for a molecule.

        Raises
        ----------
        None

        Notes
        ----------
        None

        """
        iPosition = box.coordinates[iParticle]
        wPair = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.coordinates[jParticle]
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.ForceField.cutoff2:
                wPair += self.ForceField.getPairVirial(rij2)
                if populateForces == True:
                    box.forces[iParticle] += wPair * rij / rij2
        return wPair
