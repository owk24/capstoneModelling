import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import fsolve

class ModelHelper:
    def __init__(self):
        self.R = 8.314

    def CreateFoldersForOutput(self, enzymeType) -> None:
        folderDirectoryCSVs = "./Output/{:s}/CSVs".format(enzymeType)
        folderDirectoryGraphs = "./Output/{:s}/Graphs".format(enzymeType)
        os.makedirs(Path(folderDirectoryCSVs), exist_ok=True)
        os.makedirs(Path(folderDirectoryGraphs), exist_ok=True)

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

    def GetConversionVsResidenceTimeWithCstr(self, initialS, enzymeConcentrations, Kcat, Km,
                                             listOfReactants, listOfProducts):
        def CSTR_Model(S):

            return initialS - S - residenceTime * (vMax * S) / (Km + S)

        residenceTimeArr = np.linspace(0, 480, 500)  # min
        prodFormedList = []
        substrateOutList = []
        conversionList = []

        masterDF = pd.DataFrame()
        masterDF['Residence Time'] = residenceTimeArr

        for j in enzymeConcentrations:
            vMax = Kcat * j  # (M/min)`
            for i in range(0, len(residenceTimeArr)):
                residenceTime = residenceTimeArr[i]
                finalSubstrateConcentration = fsolve(CSTR_Model, 0.03)
                if finalSubstrateConcentration < 0:
                    prodFormedList.append(prodFormedList[i - 1])
                    substrateOutList.append(substrateOutList[i - 1])
                    conversionList.append(conversionList[i - 1])
                else:
                    substrateOutList.append(finalSubstrateConcentration[0])
                    productFormed = initialS - finalSubstrateConcentration
                    prodFormedList.append(productFormed[0])
                    conversion = (initialS - finalSubstrateConcentration) / initialS
                    conversionList.append(conversion[0])

            for reactant in listOfReactants:
                masterDF["{:s} [M] - [Eo] = {:.2e}M".format(reactant, j)] = substrateOutList

            for product in listOfProducts:
                masterDF["{:s} [M] - [Eo] = {:.2e}M".format(product, j)] = prodFormedList

            masterDF["Conversion - [Eo] = {:.2e}M".format(j)] = conversionList
            prodFormedList = []
            substrateOutList = []
            conversionList = []

        return masterDF

    def GetConversionVsResidenceTimeWithCstrForCatalase(self, initialS, enzymeConcentrations, Kcat, Km):

        listOfReactants = ["H202"]
        listOfProducts = ["Water", "Oxygen"]
        def CSTR_Model(S):

            return initialS - S - residenceTime * (vMax * S) / (Km + S)
        def CSTR_Model_For_Oxygen(S):

            return initialS - S - residenceTime * (vMax * S) / (Km + S) / 2 #Becayse there is not a 1:1 ratio

        residenceTimeArr = np.linspace(0, 480, 500)  # min
        prodFormedList = []
        substrateOutList = []
        conversionList = []
        o2FormedList = []
        o2ConversionList = []

        masterDF = pd.DataFrame()
        masterDF['Residence Time'] = residenceTimeArr

        for j in enzymeConcentrations:
            vMax = Kcat * j  # (M/min)`
            for i in range(0, len(residenceTimeArr)):
                residenceTime = residenceTimeArr[i]

                finalSubstrateConcentration = fsolve(CSTR_Model, 0.03)
                if finalSubstrateConcentration < 0:
                    prodFormedList.append(prodFormedList[i - 1])
                    substrateOutList.append(substrateOutList[i - 1])
                    conversionList.append(conversionList[i - 1])
                else:
                    substrateOutList.append(finalSubstrateConcentration[0])
                    productFormed = initialS - finalSubstrateConcentration
                    prodFormedList.append(productFormed[0])
                    conversion = (initialS - finalSubstrateConcentration) / initialS
                    conversionList.append(conversion[0])

                finalSubStrateConsumptionForO2 = fsolve(CSTR_Model_For_Oxygen, 0.03)
                if finalSubStrateConsumptionForO2 < 0:
                    o2FormedList.append(o2FormedList[i - 1])
                    o2ConversionList.append(o2ConversionList[i - 1])
                else:
                    o2Formed = initialS - finalSubStrateConsumptionForO2
                    o2FormedList.append(o2Formed[0])
                    conversionO2 = (initialS - finalSubStrateConsumptionForO2) / initialS
                    o2ConversionList.append(conversionO2[0])

            masterDF["H2O2 [M] - [Eo] = {:.2e}M".format(j)] = substrateOutList
            masterDF["H2O [M] - [Eo] = {:.2e}M".format(j)] = prodFormedList
            masterDF["O2 [M] - [Eo] = {:.2e}M".format(j)] = o2FormedList
            masterDF["Conversion - [Eo] = {:.2e}M".format(j)] = conversionList
            masterDF["Conversion For O2 - [Eo] = {:.2e}M".format(j)] = o2ConversionList

            prodFormedList = []
            substrateOutList = []
            conversionList = []
            o2FormedList = []
            o2ConversionList = []

        return masterDF

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
