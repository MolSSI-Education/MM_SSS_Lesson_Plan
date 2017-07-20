#from .ForceField import LennardJones
#from .Simulation import Simulation
#from .ForceFieldManager import ForceFieldManager
#from .BoxManager import BoxManager, Box
import numpy as np
import mm_python as mmpy 

reducedTemperature = 0.851
reducedDensity = 0.776
numParticles = 100

boxLength = np.power(numParticles/reducedDensity,0.333)

print 'Box Length' , boxLength

myBox = mmpy.Box(length=boxLength)
myBoxManager = mmpy.BoxManager(myBox)

#myBoxManager.getConfigFromFile(restartFile = "initial.md.xyz", mass = 39.0)
myBoxManager.addParticles(n = numParticles, method = "lattice", mass = 39.0)

myBoxManager.assignVelocities(reducedTemperature)

myForceField = mmpy.LennardJones(cutoff= 0.5 * boxLength)

myForceFieldManager = mmpy.ForceFieldManager(myForceField)

myIntegrator = mmpy.VelocityVerlet(timeStep = 0.001, box = myBox)

mySimulation = mmpy.Simulation(
    method="molecularDynamics",
    temperature=reducedTemperature,
    steps=10000,
    printProp=1,
    printXYZ=10,
    ffManager=myForceFieldManager,
    boxManager=myBoxManager, 
    integrator = myIntegrator,
    scaleFreq = 10) 

mySimulation.run()
###
####mySimulation.getRDF(trajectory = "trajectory.xyz", bins=50)
