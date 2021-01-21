# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 21:32:35 2021

@author: Lakshay Verma and Owen Kim

"""

'''
PACKAGES
Import all relevant packages for analytics
'''

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fsolve
from Helper import modelHelper

constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")
kmOne = constDF["Km (M)"][3]
tempOne = constDF["Temperature (K)"][3]
kmTwo = constDF["Km (M)"][4]
tempTwo = constDF["Temperature (K)"][4]

invertaseHelper = modelHelper.modelHelper(kmOne, tempOne, kmTwo, tempTwo)

# activationEnergy = invertaseHelper.GetActivationEnergyFromClausiusClapeyron()
# reactionTemp = 30+273.15
# kmAtReactionTemp = invertaseHelper.GetDesiredKmFromActivationEnergy(
#     kmKnown=kmOne, tempKnown=tempOne, tempDesired=reactionTemp, activationEnergy=activationEnergy)

def substrateModelCSTR(S):
    return initialSubstrateConc - S - residenceTime*(vMax*S)/(kM+S)
initialSubstrateConc = 10000/(270*1000) #100g/L, molar mass: 270 2Da
kM = constDF["Km (M)"][3]
vMax = constDF["Vm (M/min)"][3]
residenceTime = 10 # minutes

finalSubstrateConcentration = fsolve(substrateModelCSTR, 0.005)
print(finalSubstrateConcentration)

