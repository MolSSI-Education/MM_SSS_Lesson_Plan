import numpy as np
from .particle import particle

class box(object):
    def __init__(self, length):
        self.length = length
        self.particle = []
        self.numParticles = 0

class boxManager(object):
    def __init__(self, box):
        self.box = box

    def addParticles(self, n, method):
        if method == "random":
            for iParticle in range(0, n):
                self.box.particle.append(particle())
                newPosition = (0.5 - np.random.rand(3)) * self.box.length
                self.box.particle[iParticle].position = newPosition
                self.box.numParticles = self.box.numParticles + 1

        elif method == "lattice":
            nSide = 1
            while np.power(nSide, 3) < n:
                nSide += 1
            counterPosition = np.zeros((3, ))
            for iParticle in range(0, n):
                self.box.numParticles = self.box.numParticles + 1
                self.box.particle.append(particle())
                self.box.particle[iParticle].position = \
                    (counterPosition + 0.5)*self.box.length/nSide
                counterPosition[0] += 1
                if counterPosition[0] == nSide:
                    counterPosition[0] = 0
                    counterPosition[1] += 1
                    if counterPosition[1] == nSide:
                        counterPosition[1] = 0
                        counterPosition[2] += 1

            for iParticle in range(0, self.box.numParticles):
                self.box.particle[iParticle].position -= 0.5

    def getConfigFromFile(self, restartFile):

        with open(restartFile, "r") as f:
            for lineNbr, line in enumerate(f):
                lineSplit = line.split()
                lenLine = len(lineSplit)
                if lenLine == 1:
                    self.box.numParticles = int(lineSplit[0])
                if lenLine > 1:
                    self.box.particle.append(particle())
                    iParticle = int(lineSplit[0]) - 1
                    x = float(lineSplit[1])
                    y = float(lineSplit[2])
                    z = float(lineSplit[3])

                    self.box.particle[iParticle].position[0] = x
                    self.box.particle[iParticle].position[1] = y
                    self.box.particle[iParticle].position[2] = z

                    self.box.particle[iParticle].position *= \
                            3.73
