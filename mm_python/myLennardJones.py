from .ForceField import lennardJones
from .simulation import simulation
from .ForceFieldManager import ForceFieldManager
from .boxManager import boxManager, box
import numpy as np

reducedTemperature = 0.85
reducedDensity = 0.003
numParticles = 500

sigma = 3.73
epsilon = 148.0

realTemperature = epsilon * reducedTemperature

boxLength = reducedDensity / np.power(sigma, 3)
boxLength = np.power(numParticles / boxLength, 0.333)

myBox = box(length=boxLength)
myBoxManager = boxManager(myBox)
myBoxManager.addParticles(n=numParticles, method="lattice")
#myBoxManager.getConfigFromFile(restartFile = "nistConfig.xyz")

myForceField = lennardJones(parms=(sigma, epsilon), cutoff=3 * sigma)
myForceFieldManager = ForceFieldManager(myForceField)
myForceFieldManager.assignSystemForceField(myBox)

mySimulation = simulation(
    method="monteCarlo",
    temperature=realTemperature,
    steps=1000000,
    printProp=1000,
    printXYZ=1000,
    maxDisp=sigma,
    ffManager=myForceFieldManager,
    boxManager=myBoxManager)

mySimulation.run()

#mySimulation.analyze.getRDF(parms)
