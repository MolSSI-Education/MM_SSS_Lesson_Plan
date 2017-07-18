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

    def getTotalPairEnergy(self, box):
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
        ePair = 0.0
        for iParticle in range(0, box.numParticles):
            eInter = self.getMolEnergy(iParticle, box)
            ePair = ePair + eInter
        ePair = ePair/2.0
        return ePair

    def getMolVirial(self, iParticle, box):
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
        return wPair

    def getSystemVirial(self, box):
        """
	Computes the system virial.

        Parameters
        ----------
	box: box
	The box containing the particles.

        Returns
        ----------
        wTotal: float
	Total intermolecular virial.

        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        wTotal = 0.0
        for iParticle in range(0, box.numParticles):
            wInter = self.getMolVirial(iParticle, box)
            wTotal = wTotal + wInter
        wTotal = wTotal/2.0
        return wTotal
