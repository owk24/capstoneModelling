import numpy as np

class modellingHelper():
    def __init__(self, kmOne, tempOne, kmTwo, tempTwo):
        self.R = 8.314
        self.kmOne = kmOne
        self.tempOne = tempOne
        self.kmTwo = kmTwo
        self.tempTwo = tempTwo

    def ClausiusClapeyron(self) -> float:
        activationEnergy = self.R*(self.tempOne*self.tempTwo/(self.tempTwo - self.tempOne))*np.log(self.kmOne/self.kmTwo)

        return activationEnergy

    def BindingAffinity(self, activationEnergy) -> float:
        kmTwo = self.kmOne*np.exp((-activationEnergy/self.R)((1 / self.tempTwo) - (1 / self.tempOne)))

        return kmTwo

    def ModelReaction(self, vMax, S, Km) -> float:
        reactionRate = vMax * (S / (Km + S))

        return reactionRate
