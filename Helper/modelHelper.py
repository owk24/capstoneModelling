import numpy as np

class modelHelper():
    def __init__(self, kmOne, tempOne, kmTwo, tempTwo):
        self.R = 8.314
        self.kmOne = kmOne
        self.tempOne = tempOne
        self.kmTwo = kmTwo
        self.tempTwo = tempTwo

    def ClausiusClapeyronReturnActivationEnergy(self) -> float:
        activationEnergy = self.R * (self.tempOne*self.tempTwo/(self.tempTwo - self.tempOne)) \
                           * np.log(self.kmOne / self.kmTwo)

        return activationEnergy

    def BindingAffinityReturnDesiredKm(self, kmKnown, tempKnown, tempDesired, activationEnergy) -> float:
        kmDesired = kmKnown * np.exp((-activationEnergy/self.R) * ((1/tempDesired) - (1/tempKnown)))

        return kmDesired

    def MichaelisModelReactionReturnRate(self, vMax, S, Km) -> float:
        reactionRate = vMax * S / (Km + S)
        
        return reactionRate
