import numpy as np


class ForceFieldManager(object):
    def __init__(self, ForceField):
        self.ForceField = ForceField

    def getMolPairEnergyAndVirial(self, iParticle, box, populateForces = False):

        iPosition = box.coordinates[iParticle]
        eTotal = 0.0
        wTotal = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.coordinates[jParticle]
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.ForceField.cutoff2:
                ePair = self.ForceField.evaluate(rij2)
                wPair = self.ForceField.getPairVirial(rij2)

                eTotal += ePair
                wTotal += wPair
                if populateForces == True:
                    box.forces[iParticle] += wPair * (rij / rij2)

        return eTotal, wTotal

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
