from forceField import lennardJones
from box import box
from simulation import simulation
from forceFieldManager import forceFieldManager
from boxManager import boxManager
import numpy as np

myBox = box(length = 37.3)
myBoxManager = boxManager(myBox)
#myBoxManager.addParticles(n=200,method="lattice") 
myBoxManager.getConfigFromFile(restartFile = "lj_sample_config_periodic1.txt")

myForceField = lennardJones(parms = (3.73, 148.0), cutoff = 11.19)
myForceFieldManager = forceFieldManager(myForceField)
myForceFieldManager.assignSystemForceField(myBox)

mySimulation = simulation(method = "monteCarlo", temperature = 133.2, 
        steps = 1000000, printProp = 1000, printXYZ = 1000, maxDisp = 0.1,
        ffManager = myForceFieldManager, boxManager = myBoxManager)

mySimulation.run()

#mySimulation.analyze.getRDF(parms)
