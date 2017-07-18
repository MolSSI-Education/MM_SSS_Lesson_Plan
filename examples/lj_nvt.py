#from .ForceField import LennardJones
#from .Simulation import Simulation
#from .ForceFieldManager import ForceFieldManager
#from .BoxManager import BoxManager, Box
import numpy as np
import mm_python as mmpy 

reducedTemperature = 0.9
reducedDensity = 0.9
numParticles = 100

sigma = 3.73
epsilon = 148.0

realTemperature = epsilon * reducedTemperature

boxLength = reducedDensity / np.power(sigma, 3)
boxLength = np.power(numParticles / boxLength, 0.333)

print 'Box Length' , boxLength

myBox = mmpy.Box(length=boxLength)
myBoxManager = mmpy.BoxManager(myBox)
myBoxManager.addParticles(n=numParticles, method="lattice")
#myBoxManager.getConfigFromFile(restartFile = "initial.xyz")

myForceField = mmpy.LennardJones(
    parms=(sigma, epsilon), cutoff= 0.5 * boxLength)

myForceFieldManager = mmpy.ForceFieldManager(myForceField)

mySimulation = mmpy.Simulation(
    method="monteCarlo",
    temperature=realTemperature,
    steps=50000,
    printProp=1000,
    printXYZ=1000,
    maxDisp=0.1*sigma,
    ffManager=myForceFieldManager,
    boxManager=myBoxManager)

mySimulation.run()

#mySimulation.getRDF(trajectory = "trajectory.xyz", bins=50)
