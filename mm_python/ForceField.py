import numpy as np

class ForceField(object):
    """
    Base class for Force Fields. The attributes are the
    type of the potential (e.g. LJ, Mie), energy cut off (cutoff)
    and the energy cut off squared (cutoff2)

    Properties
    ----------
    cutoff : float
        Cutoff distance (Angstroms)
    cutoff2 : float
        Squared cutoff (Angstroms^2)
    switch : float or None
        If not None, the distance at which the switch is to be turned on.

    Parameters
    ----------
    potentialType : str
        String specifying the potential energy type (e.g. "lennardJones")
    cutoff : float
        Custoff distance (Anstroms)
    switch : float, optional, default=None
        If not `None`, the distance at which the switch will be turned on.
        If `None`, no switch is used; only a cutoff.

    """
    def __init__(self, potentialType, cutoff, switch=None):

        self.potentalType = potentialType
        self.cutoff = cutoff
        self.cutoff2 = np.power(cutoff, 2)
        if switch:
            if (switch > cutoff):
                raise Exception("'switch' must be in range [0,cutoff]")

            self.switch = switch
            self.switch2 = np.power(switch, 2)
        else:
            self.switch = None
            self.switch2 = None

class LennardJones(ForceField):
    """
    Class for simple pairwise Lennard Jones potentials. This class is a
    child of the base class ForceField. Two functions are defined: evaluate
    getPairVirial, and getTailCorrection and getPressureCorrection.
    """

    def __init__(self, cutoff, switch=None):
        """
        cutoff : float
            Custoff distance (Anstroms)
        switch : float, optional, default=None
            If not `None`, the distance at which the switch will be turned on.
            If `None`, no switch is used; only a cutoff.
        """
        ForceField.__init__(self, "lennardJones", cutoff, switch)

    def evaluate(self, rij2):
        """
        Computes the interaction energy between two particles using
        the Lennard Jones potential.

        Parameters
        ----------
        rij2: float
            Squared distance between two particles (Anstroms^2)

        Returns
        ----------
        energy: float
            Interaction potential energy between two particles (kJ/mol?)

        Raises
        ----------
        None

        Notes
        ----------
        None

        """
        sigByR2 = 1.0 / rij2
        sigByR6 = np.power(sigByR2, 3)
        sigByR12 = np.power(sigByR6, 2)
        return 4.0 * (sigByR12 - sigByR6)

    def computeSwitch(self, rij2):
        """
        Compute the value of the switch function S(r).

        S(r) = 1 if r <= switch, 0 if r >= cutoff, otherwise
        S(r) = (r_cut^2 - r_ij^2)^2 (r_cut^2 + 2 r_ij^2 - 3 r_sw^2) / (r_cut^2 - r_sw^2)^3

        Parameters
        ----------
        rij2 : float
            Squared distance between atoms i and j (Angstroms^2)

        Returns
        -------
        S : float
            Switch function S(r)

        """
        if self.switch is None:
            # No switch
            if rij2 < self.cutoff2:
                S = 1.0
            else:
                S = 0.0
        else:
            # Switch in use
            if rij2 <= self.switch2:
                S = 1.0
            elif rij2 >= self.cutoff2:
                S = 0.0
            else:
                S = (self.cutoff2 - rij2)**2 * (self.cutoff2 + 2*rij2 - 3*self.switch2) / (self.cutoff2 - self.switch2)**3

        return S

    def computeSwitchDerivative(self, rij2):
        """
        Compute the derivative of the switch function, dS/dr.

        Parameters
        ----------
        rij2 : float
            Squared distance between atoms i and j (Angstroms^2)

        Returns
        -------
        dS_dr : float
            Scalar derivative dS/dr

        """
        if self.switch is None:
            # No switch
            dS_dr = 0.0
        else:
            # Switch in use
            if (rij2 <= self.switch2) or (rij2 >= self.cutoff2):
                dS_dr = 0.0
            else:
                dS_dr = 12 * np.sqrt(rij2) * (rij2 - self.cutoff2) * (rij2 - self.switch2) / (self.cutoff2 - self.switch2)**3

        return dS_dr

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

        sigByR2 =  1 / rij2
        sigByR6 = np.power(sigByR2,3)
        sigByR12 = np.power(sigByR6,2)
        pairVirial = 24.0 * (2.0*sigByR12 - sigByR6)
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
        sigByCutoff3 = np.power(1.0/self.cutoff, 3)
        sigByCutoff9 = np.power(sigByCutoff3, 3)
        eCorrection = sigByCutoff9 - 3.0 * sigByCutoff3
        eCorrection = 8.0/9.0 * np.pi \
                * box.numParticles/np.power(box.length, 3) \
                * box.numParticles \
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
                np.power(1.0/self.cutoff, 3)
        sigByCutoff9 = np.power(sigByCutoff3, 3)
        pCorrection = 2.0/3.0*sigByCutoff9 - sigByCutoff3
        pCorrection = 16.0/3.0 * np.pi \
                * np.power(box.numParticles/np.power(box.length, 3),2) \
                * pCorrection
        return pCorrection
