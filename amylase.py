# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 23:32:20 2020

@author: Lakshay Verma and Owen Kim

"""

'''
PACKAGES
Import all relevant packages for analytics
'''

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import fsolve
from Helper.ModelHelper import ModelHelper

def substrateModelCSTR(S):
    return initSubstrateConc - S - residenceTime*(vMax*S)/(bindingAffinity+S)

constDF = pd.read_csv(Path("Constants") / "amylaseConstants.csv")   
dfToArr = constDF.to_numpy()
molarMassOfProduct = 180.16 #g/mol

amylaseHelper = ModelHelper()

vectorModel = np.vectorize(amylaseHelper.GetRateFromMichaelisModelReaction)
amylaseReactions = []

S_amylase = np.linspace(2e-7, 6e-6, 100)

for i in range(5, 10):
    amylaseReactions.append(vectorModel(dfToArr[i][1], S_amylase, 
                                        dfToArr[i][0]))

fig = plt.figure()
for i in range(5, 10):
    plt.plot(S_amylase, amylaseReactions[i-5], label=str(dfToArr[i][4]) + 
             '($^\circ$C)')

plt.xlabel('Substrate Concentration [M]')
plt.ylabel('Reaction Rate [M/min]')
plt.title('Modelling of Amylase Reaction', fontweight='bold')
plt.legend()
plt.show()

#Use Clausius Clapeyron with given datapoints to determine activation energy
#NOTE: This is not in use anymore. Km and Vmax were found at optimal conditions
#via a paper, and thus no need to extrapolate values.
activationEnergy = amylaseHelper.GetActivationEnergyFromClausiusClapeyron(
                                        tempOne=dfToArr[5][5],
                                        kmOne=dfToArr[5][0],
                                        tempTwo=dfToArr[9][5],
                                        kmTwo=dfToArr[9][0])
reactorTemp = 70+273.15

kmAtReactorConditions = amylaseHelper.GetDesiredKmFromActivationEnergy(
    kmKnown=dfToArr[5][0], tempKnown=dfToArr[5][5],
    tempDesired = reactorTemp, activationEnergy=activationEnergy)

'''
MASS BALANCE
'''

initialMassStarch = 5 # g Starch/L
initSubstrateConc = initialMassStarch/359.33 # mol Starch/L 
initEnzymeConcentration = 9e-8 #M
kcat = constDF['Kcat (1/min)'][9] #1/min
vMax = kcat * initEnzymeConcentration #(M/min)
bindingAffinity = constDF['Km (M)'][9] #M

residenceTime = 120

finalSubstrateConcentration = fsolve(substrateModelCSTR, 0.03)
productFormed = initSubstrateConc - finalSubstrateConcentration
gramsOfProductFormed = productFormed * molarMassOfProduct
print("\nConditions: " + str(constDF['Temperature (C)'][9]) + " degrees Celsius")
print("For a residence time of " + str(residenceTime) + " minutes, and an initial" 
      + " substrate concentration of " + str(initialMassStarch) + " g starch/L, " +
      str(gramsOfProductFormed) + " g/L of glucose formed")
