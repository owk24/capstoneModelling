# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 23:32:20 2020

@author: Lakshay Verma and Owen Kim

"""

'''
PACKAGES
Import all relevant packages for analytics
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import fsolve
from Helper import modelHelper

def substrateModelCSTR(S):
    return initSubstrateConc - S - residenceTime*(vMax*S)/(bindingAffinity+S)

amylaseHelper = modelHelper.modelHelper(kmOne=69, tempOne=69, kmTwo=69, tempTwo=69)

constDF = pd.read_csv('constantValues.csv')    
dfToArr = constDF.to_numpy()
vectorModel = np.vectorize(amylaseHelper.GetRateFromMichaelisModelReaction)
amylaseReactions = []

S_amylase = np.linspace(0.0001, 0.00005, 100)

for i in range(0, 5):
    amylaseReactions.append(vectorModel(dfToArr[4][i+2], S_amylase, dfToArr[2][i+2]))

fig = plt.figure()
for i in range(0, 5):
    plt.plot(S_amylase, amylaseReactions[i], label=str(dfToArr[0][i+2]))

plt.xlabel('Substrate Concentration [M]')
plt.ylabel('Reaction Rate [M/min]')
plt.title('Modelling of Amylase Reaction', fontweight='bold')
plt.legend()
plt.show()

'''
MASS BALANCE
'''
 
initSubstrateConc = 5/58400
flowrate = 5 #(L/min)
volumeOfReactor = 5 #L
vMax = 2.3004599999999998e-08 #(M/min)
bindingAffinity = 4.29e-05 #(M)

residenceTime = volumeOfReactor/flowrate

substrateModellingEquation = fsolve(substrateModelCSTR, 0.000001)
