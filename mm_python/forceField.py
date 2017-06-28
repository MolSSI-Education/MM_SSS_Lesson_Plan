"""
This is a file for my forcefields
"""

import numpy as np


class ForceField(object):
    def __init__(self, potentialType, parms, cutoff):
        """
        This is a base class for my force fields.
        """

        if not isinstance(parms, tuple):
            raise TypeError("ForceField: second parameters, params, must be passed as a tuple.")

        self.parms = parms
        self.potentalType = potentialType
        self.cutoff = cutoff
        self.cutoff2 = np.power(cutoff, 2)


class LennardJones(ForceField):
    def __init__(self, parms, cutoff):
        ForceField.__init__(self, "lennardJones", parms, cutoff)

    def evaluate(self, rij2):
        sigByR2 = np.power(self.parms[0], 2) / rij2
        sigByR6 = np.power(sigByR2, 3)
        sigByR12 = np.power(sigByR6, 2)
        return 4.0 * self.parms[1] * (sigByR12 - sigByR6)

    def getTailCorrection(self, box):
        sigByCutoff3 = \
                np.power(self.parms[0]/self.cutoff,3)
        sigByCutoff9 = np.power(sigByCutoff3, 3)
        eCorrection = sigByCutoff9 - 3.0 * sigByCutoff3
        eCorrection = 8.0/9.0 * np.pi * \
                box.numParticles/np.power(box.length,3) * \
                np.power(self.parms[0],3) * \
                box.numParticles * \
                self.parms[1] * \
                eCorrection
        return eCorrection
