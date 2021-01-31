import numpy as np
import pandas as pd
from scipy.integrate import odeint

class ModelHelper:
    def __init__(self):
        self.R = 8.314

    def GetActivationEnergyFromClausiusClapeyron(self, tempOne, kmOne, tempTwo, kmTwo) -> float:
        activationEnergy = self.R * (self.tempOne*self.tempTwo/(self.tempTwo - self.tempOne)) \
                           * np.log(self.kmOne / self.kmTwo)

        return activationEnergy

    def GetDesiredKmFromActivationEnergy(self, kmKnown, tempKnown, tempDesired, activationEnergy) -> float:
        kmDesired = kmKnown * np.exp((-activationEnergy/self.R) * ((1/tempDesired) - (1/tempKnown)))

        return kmDesired

    def GetRateFromMichaelisModelReaction(self, vMax, S, Km) -> float:
        reactionRate = vMax * S / (Km + S)
        
        return reactionRate

    def GetRatesVersusTimeFromModel(self, t, km, kcat, initialConditions, initialEnzymeConc):
        def model(inputs, t, Km, Kcat, initialEnzymeConcentration):
            Vm = Kcat * initialEnzymeConcentration
            r = Vm * inputs[2] / (Km + inputs[2])
            dGdt = r
            dFdt = r
            dSdt = -r
            outputs = [dGdt, dFdt, dSdt]

            return outputs
        concentrationOfComponents = odeint(model, initialConditions, t, args=(km, kcat, initialEnzymeConc))
        concentrationDf = pd.DataFrame({"Glucose [M/min]": concentrationOfComponents[:, 0],
                                        "Fructose [M/min]": concentrationOfComponents[:, 1],
                                        "Sucrose [M/min]": concentrationOfComponents[:, 2]})
        concentrationDf.index.name = "Time [min]"

        return concentrationDf
