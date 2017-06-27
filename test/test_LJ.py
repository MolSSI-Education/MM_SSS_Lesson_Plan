
import numpy as np
import pytest
import mm_python as mmpy


def test_lj():

    params = (2.5, 6.5)
    cutoff = 10.0

    lj_ff = mmpy.LennardJones(params, cutoff) 

    rij2 = np.linspace(0.1, 50, 1000)
    
    bench_energy =  (np.power(params[0], 2) / rij2) ** 6 
    bench_energy -= (np.power(params[0], 2) / rij2) ** 3
    bench_energy *= 4.0 * params[1]

    energy = lj_ff.evaluate(rij2)

    assert np.allclose(energy, bench_energy)
