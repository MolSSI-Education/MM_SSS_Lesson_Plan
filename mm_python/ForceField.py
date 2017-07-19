import numpy as np

class ForceField(object):
    """
    Base class for Force Fields. The main attributes are the
    potential parameters (parms)
    type of the potential (e.g. LJ, Mie), energy cut off (cutoff)
    and the energy cut off squared (cutoff2)
    """
    def __init__(self, potentialType, parms, cutoff):

        if not isinstance(parms, tuple):
            raise TypeError("ForceField: second parameters, params, \
				must be passed as a tuple.")

        self.parms = parms
        self.potentalType = potentialType
        self.cutoff = cutoff
        self.cutoff2 = np.power(cutoff, 2)


class LennardJones(ForceField):
    """
    Class for simple pairwise Lennard Jones potentials. This class is a 
    child of the base class ForceField. Two functions are defined: evaluate
    and getTailCorrections. 
    """

    def __init__(self, parms, cutoff):
        ForceField.__init__(self, "lennardJones", parms, cutoff)

    def evaluate(self, rij2):
        """
        Computes the interaction energy between two particles using
	the Lennard Jones potential.
            
        Parameters
        ----------
	rij2: float
	Distance between two particles

        Returns
        ----------
        energy: float
	Interaction energy between two particles

        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        sigByR2 = np.power(self.parms[0], 2) / rij2
        sigByR6 = np.power(sigByR2, 3)
        sigByR12 = np.power(sigByR6, 2)
        return 4.0 * self.parms[1] * (sigByR12 - sigByR6)


    def getPairVirial(self, rij2):
        """
        Computes the virial between particles i and j.
            
        Parameters
        ----------
	rij2: float
	Distance between two particles

        Returns
        ----------
        pairVirial: float
        Force between two particles

        Raises
        ----------
        None


        Notes
        ----------
        None

        """

	sigByR2 =  np.power(self.parms[0],2) / rij2
        sigByR6 = np.power(sigByR2,3)
        sigByR12 = np.power(sigByR6,2)
        pairVirial = 24.0 * self.parms[1] * (2.0*sigByR12 - sigByR6)
	return pairVirial


    def getTailCorrection(self, box):
	"""
	Computes the system tail correction due to energy cut off.            
 
        Parameters
        ----------
	box: box
	Box for which the tail correction will be computed. 

        Returns
        ----------
	eCorrection: float
	The computed tail correction.

        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        sigByCutoff3 = \
                np.power(self.parms[0]/self.cutoff, 3)
        sigByCutoff9 = np.power(sigByCutoff3, 3)
        eCorrection = sigByCutoff9 - 3.0 * sigByCutoff3
        eCorrection = 8.0/9.0 * np.pi \
                * box.numParticles/np.power(box.length, 3) \
                * np.power(self.parms[0], 3) \
                * box.numParticles \
                * self.parms[1] \
                * eCorrection
        return eCorrection

    def getPressureCorrection(self, box):

	"""
	Computes the pressure correction due to cut off.            
 
        Parameters
        ----------
	box: box
	Box for which the tail correction will be computed. 

        Returns
        ----------
	pCorrection: float
	The computed pressure tail correction.

        Raises
        ----------
        None


        Notes
        ----------
        None

        """
        sigByCutoff3 = \
                np.power(self.parms[0]/self.cutoff, 3)
        sigByCutoff9 = np.power(sigByCutoff3, 3)
        pCorrection = 2.0/3.0*sigByCutoff9 - sigByCutoff3
        pCorrection = 16.0/3.0 * np.pi \
                * np.power(box.numParticles/np.power(box.length, 3),2) \
                * np.power(self.parms[0], 3) \
                * self.parms[1] \
                * pCorrection
        return pCorrection
