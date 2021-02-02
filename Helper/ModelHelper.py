import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import odeint

class ModelHelper:
    def __init__(self):
        self.R = 8.314

    def CreateFoldersForOutput(self, enzymeType) -> None:
        folderDirectoryCSVs = "./{:s}/CSVs".format(enzymeType)
        folderDirectoryGraphs = "./{:s}/Graphs".format(enzymeType)
        os.makedirs(Path(folderDirectoryCSVs), exist_ok=True)
        os.makedirs(Path(folderDirectoryGraphs).format(enzymeType), exist_ok=True)

    def GetActivationEnergyFromClausiusClapeyron(self, tempOne, kmOne, tempTwo, kmTwo) -> float:
        activationEnergy = self.R * (tempOne*tempTwo/(tempTwo - tempOne)) \
                           * np.log(kmOne / kmTwo)

        return activationEnergy

    def GetDesiredKmFromActivationEnergy(self, kmKnown, tempKnown, tempDesired, activationEnergy) -> float:
        kmDesired = kmKnown * np.exp((-activationEnergy/self.R) * ((1/tempDesired) - (1/tempKnown)))

        return kmDesired

    def GetRateFromMichaelisModelReaction(self, Km, Vmax, S) -> float:
        reactionRate = Vmax * S / (Km + S)
        
        return reactionRate

    def GetRatesVersusTimeFromModel(self, t, Km, Vm, initialConditions):
        def model(inputs, t, Km, Vm):
            r = Vm * inputs[2] / (Km + inputs[2])
            dGdt = r
            dFdt = r
            dSdt = -r
            outputs = [dGdt, dFdt, dSdt]

            return outputs

        concentrationOfComponents = odeint(model, initialConditions, t, args=(Km, Vm))
        forwardRateArray = self.GetRateFromMichaelisModelReaction(Km, Vm, concentrationOfComponents[:, 2])
        concentrationDf = pd.DataFrame({"Time [min]": t,
                                        "Glucose [M]": concentrationOfComponents[:, 0],
                                        "Fructose [M]": concentrationOfComponents[:, 1],
                                        "Sucrose [M]": concentrationOfComponents[:, 2],
                                        "Forward Rate [M/s]": forwardRateArray})

        return concentrationDf
