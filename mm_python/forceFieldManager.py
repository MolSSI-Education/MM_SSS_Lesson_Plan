import numpy as np


class forceFieldManager(object):
    def __init__(self, forceField):
        self.forceField = forceField

    def assignSystemForceField(self, box):
        for iParticle in range(0, box.numParticles):
            box.particle[iParticle].parms = self.forceField.parms

    def getMolEnergy(self, iParticle, box):
        iPosition = box.particle[iParticle].position
        eij = 0.0
        for jParticle in range(0, box.numParticles):
            if iParticle == jParticle: continue
            jPosition = box.particle[jParticle].position
            rij = iPosition - jPosition
            rij = rij - box.length * np.round(rij / box.length)
            rij2 = np.sum(np.power(rij, 2))
            if rij2 < self.forceField.cutoff2:
                eij += self.forceField.evaluate(rij2)
        return eij

    def getPairEnergy(self, box):
        ePair = 0.0
        for iParticle in range(0, box.numParticles):
            eInter = self.getMolEnergy(iParticle, box)
            box.particle[iParticle].energy = eInter
            ePair += eInter / 2.0
        ePair = ePair
        return ePair
