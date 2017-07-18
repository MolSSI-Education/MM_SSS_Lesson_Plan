
import numpy as np
import pytest
import mm_python as mmpy
import matplotlib.pyplot as plt

def test_force():

    params = (3.73,148.0)
    cutoff = 14.0
    
    lj_ff = mmpy.LennardJones(params, cutoff) 
    
    rij = np.linspace(3.5, 20 , 100)
    
    bench_force =  2.0*np.power(params[0] / rij,12)
    bench_force -= np.power(params[0] / rij,6)
    bench_force *= 24.0 * params[1] / rij
    
    rij2 = np.power(rij,2)
    force = lj_ff.getPairVirial(rij2)/rij
    
    assert np.allclose(force, bench_force)
