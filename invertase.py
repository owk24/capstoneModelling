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
from scipy.integrate import odeint
from Helper import modelHelper

constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")
kmOne = constDF["Km (M)"][3]
tempOne = constDF["Temperature (K)"][3]
kmTwo = constDF["Km (M)"][4]
tempTwo = constDF["Temperature (K)"][4]

invertaseHelper = modelHelper.modelHelper(kmOne, tempOne, kmTwo, tempTwo)

def model(inputs, t):
    index = 3
    Vm = constDF["Vm (M/min)"][index]
    Km = constDF["Km (M)"][index]
    r = Vm*inputs[2]/(Km + inputs[2])
    dGdt = r
    dFdt = r
    dSdt = -r
    outputs = [dGdt, dFdt, dSdt]
    return outputs

initialG = 0
initialF = 0
initialS = 100/58400 #100g/L, 58400kDa molar mass

initialConditions = [initialG, initialF, initialS]

t = np.linspace(0, 50, 1000)

odeModel = odeint(model, initialConditions, t)

plt.plot(t, odeModel[:, 0], label="Glucose")
plt.plot(t, odeModel[:, 1], label="Fructose")
plt.plot(t, odeModel[:, 2], label="Sucrose")
plt.xlabel("Time (min)")
plt.ylabel("Concentration (M)")
plt.legend(loc="best")
plt.show()

# def substrateModelCSTR(S):
#     return initialSubstrateConc - S - residenceTime*(vMax*S)/(kM+S)
# initialSubstrateConc = 10000/(270*1000) #100g/L, molar mass: 270 2Da
# kM = constDF["Km (M)"][3]
# vMax = constDF["Vm (M/min)"][3]
# residenceTime = 10 # minutes

# finalSubstrateConcentration = fsolve(substrateModelCSTR, 0.005)
# print(finalSubstrateConcentration)

