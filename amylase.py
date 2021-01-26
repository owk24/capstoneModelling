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
from Helper import modelHelper

def substrateModelCSTR(S):
    return initSubstrateConc - S - residenceTime*(vMax*S)/(bindingAffinity+S)

constDF = pd.read_csv(Path("Constants") / "amylaseConstants.csv")   
dfToArr = constDF.to_numpy()

amylaseHelper = modelHelper.modelHelper(kmOne=dfToArr[5][0], tempOne=dfToArr[5][5], 
                                        kmTwo=dfToArr[9][0], tempTwo=dfToArr[9][5])

vectorModel = np.vectorize(amylaseHelper.GetRateFromMichaelisModelReaction)
amylaseReactions = []

S_amylase = np.linspace(2e-7, 5e-7, 100)

for i in range(5, 10):
    amylaseReactions.append(vectorModel(dfToArr[i][1], S_amylase, dfToArr[i][0]))

fig = plt.figure()
for i in range(5, 10):
    plt.plot(S_amylase, amylaseReactions[i-5], label=str(dfToArr[i][4]) + '($^\circ$C)')

plt.xlabel('Substrate Concentration [M]')
plt.ylabel('Reaction Rate [M/min]')
plt.title('Modelling of Amylase Reaction', fontweight='bold')
plt.legend()
plt.show()

#Use Clausius Clapeyron with given datapoints to determine activation energy

activationEnergy = amylaseHelper.GetActivationEnergyFromClausiusClapeyron()
reactorTemp = 70+273.15
kmAtReactorConditions = amylaseHelper.GetDesiredKmFromActivationEnergy(kmKnown=dfToArr[5][0], tempKnown=dfToArr[5][5],
                                                                       tempDesired = reactorTemp, activationEnergy=activationEnergy)
print(kmAtReactorConditions)

'''
MASS BALANCE
'''

initSubstrateConc = 5/359.33 # 5 g/L Starch
#flowrate = 5 #(L/min)
#volumeOfReactor = 5 #L
vMax = dfToArr[10][1] #(M/min)
bindingAffinity = dfToArr[10][0] #(M)

residenceTime = 240

substrateModellingEquation = fsolve(substrateModelCSTR, 0.000001)
