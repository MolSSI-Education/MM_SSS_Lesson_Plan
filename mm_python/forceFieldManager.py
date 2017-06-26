import numpy as np


class ForceFieldManager(object):
    def __init__(self, ForceField):
        self.ForceField = ForceField

    def assignSystemForceField(self, box):
        for iParticle in range(0, box.numParticles):
            box.particle[iParticle].parms = self.ForceField.parms

    def getMolEnergy(self, iParticle, box):
        iPosition = box.particle[iParticle].position
        eij = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.particle[jParticle].position
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.ForceField.cutoff2:
                eij += self.ForceField.evaluate(rij2)
        return eij

    def getPairEnergy(self, box):
        ePair = 0.0
        for iParticle in range(0, box.numParticles):
            eInter = self.getMolEnergy(iParticle, box)
            box.particle[iParticle].energy = eInter
            ePair += eInter / 2.0
        ePair = ePair
        return ePair
