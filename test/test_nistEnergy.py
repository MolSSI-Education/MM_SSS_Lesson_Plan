import numpy as np
import pytest
import mm_python as mmpy


def test_nistEnergy():


    answer = (-4351.5,-198.49)

    params = (3.73, 148.0)
    rCut = 3*params[0]

    myBox = mmpy.Box(length=10.0*params[0])
    myBoxManager = mmpy.BoxManager(myBox)
    myBoxManager.getConfigFromFile(restartFile = "test/nistConfig.xyz")

    myForceField = mmpy.LennardJones(parms = params, cutoff = rCut)
    ffManager = mmpy.ForceFieldManager(myForceField)

    totalEnergy = ffManager.getTotalEnergy(myBox)
    correctionEnergy = myForceField.getTailCorrection(myBox)

    bench_energy = (totalEnergy/params[1], correctionEnergy/params[1])
    assert np.allclose(answer, bench_energy)

