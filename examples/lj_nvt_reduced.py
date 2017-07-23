#from .ForceField import LennardJones
#from .Simulation import Simulation
#from .ForceFieldManager import ForceFieldManager
#from .BoxManager import BoxManager, Box
import numpy as np
import mm_python as mmpy

reducedTemperature = 0.9
reducedDensity = 0.9
numParticles = 100


boxLength = np.power(numParticles/reducedDensity,0.333)

print(boxLength)

myBox = mmpy.Box(length=boxLength)
myBoxManager = mmpy.BoxManager(myBox)
myBoxManager.addParticles(n=numParticles, method="lattice")

myForceField = mmpy.LennardJones(cutoff= 0.5 * boxLength)

myForceFieldManager = mmpy.ForceFieldManager(myForceField)

mySimulation = mmpy.Simulation(
    method="monteCarlo",
    temperature=reducedTemperature,
    steps=50000,
    printProp=1000,
    printXYZ=1000,
    maxDisp=0.1,
    ffManager=myForceFieldManager,
    boxManager=myBoxManager)

mySimulation.run()

#mySimulation.getRDF(trajectory = "trajectory.xyz", bins=50)
