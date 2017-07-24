import numpy as np
import pytest
import mm_python as mmpy


def test_nistEnergy():


    answer = (-4351.5,-198.49)

    rCut = 3.0

    myBox = mmpy.Box(length=10.0)
    myBoxManager = mmpy.BoxManager(myBox)
    myBoxManager.getConfigFromFile\
            (restartFile = "test/lj_sample_config_periodic1.txt")

    myForceField = mmpy.LennardJones(cutoff = rCut)
    ffManager = mmpy.ForceFieldManager(myForceField)

    totalPairEnergy, wPair = ffManager.getTotalPairEnergyAndVirial(myBox)
    correctionEnergy = myForceField.getTailCorrection(myBox)

    bench_energy = (totalPairEnergy, correctionEnergy)
    assert np.allclose(answer, bench_energy)

