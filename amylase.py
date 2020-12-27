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
from Helper import modellingHelper

amylaseHelper = modellingHelper.modellingHelper(kmOne=69, kmTwo=69)

def modelReaction(kCat, vMax, S, Km, *args):
    reactionRate = vMax*(S/(Km+S))
    return reactionRate

constDF = pd.read_csv('constantValues.csv')    
dfToArr = constDF.to_numpy()
vectorModel = np.vectorize(modelReaction)
S = np.linspace(0, 0.5, 100)

amylaseReactions = []

S_amylase = np.linspace(0.0000004, 0.0000005, 100)

for i in range(0, 5):
    amylaseReactions.append(vectorModel(dfToArr[3][i+2], dfToArr[4][i+2], S_amylase, dfToArr[2][i+2]))

fig = plt.figure()
for i in range(0, 5):
    plt.plot(S_amylase, amylaseReactions[i], label=str(dfToArr[0][i+2]) + '($^\circ$C)')


plt.xlabel('Substrate Concentration [M]')
plt.ylabel('Reaction Rate [M/min]')
plt.title('Modelling of Amylase Reaction', fontweight='bold')
plt.legend()
plt.show()

