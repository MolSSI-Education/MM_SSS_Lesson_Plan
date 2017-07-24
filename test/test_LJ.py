
import numpy as np
import pytest
import mm_python as mmpy


def test_lj():

    cutoff = 10.0

    lj_ff = mmpy.LennardJones(cutoff) 

    rij2 = np.linspace(0.1, 50, 1000)
    
    bench_energy =  (1.0 / rij2) ** 6 
    bench_energy -= (1.0 / rij2) ** 3
    bench_energy *= 4.0 

    energy = lj_ff.evaluate(rij2)

    assert np.allclose(energy, bench_energy)

def test_switch():
    """
    Test that switch function has correct limits.

    """

    switch = 9.0
    cutoff = 10.0

    lj_ff = mmpy.LennardJones(params, cutoff, switch=switch)

    assert lj_ff.computeSwitch(switch**2) == 1.0
    assert lj_ff.computeSwitch(cutoff**2) == 0.0

    assert lj_ff.computeSwitchDerivative(switch**2) == 0.0
    assert lj_ff.computeSwitchDerivative(((switch+cutoff)/2.0)**2) != 0.0
    assert lj_ff.computeSwitchDerivative(cutoff**2) == 0.0
