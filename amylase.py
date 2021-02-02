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
from Helper.ModelHelper import ModelHelper

constDF = pd.read_csv(Path("Constants") / "AmylaseConstants.csv")   
dfToArr = constDF.to_numpy()
molarMassOfProduct = 180.16 #g/mol

amylaseHelper = ModelHelper()
amylaseHelper.CreateFoldersForOutput("Amylase")

# vectorModel = np.vectorize(amylaseHelper.GetRateFromMichaelisModelReaction)
# amylaseReactions = []
#
# S_amylase = np.linspace(2e-7, 6e-6, 100)
#
# for i in range(5, 10):
#     amylaseReactions.append(vectorModel(dfToArr[i][1], S_amylase,
#                                         dfToArr[i][0]))
#
# fig = plt.figure()
# for i in range(5, 10):
#     plt.plot(S_amylase, amylaseReactions[i-5], label=str(dfToArr[i][4]) +
#              '($^\circ$C)')
#
# plt.xlabel('Substrate Concentration [M]')
# plt.ylabel('Reaction Rate [M/min]')
# plt.title('Modelling of Amylase Reaction', fontweight='bold')
# plt.legend()
# plt.show()
#
# # Use Clausius Clapeyron with given datapoints to determine activation energy
# # NOTE: This is not in use anymore. Km and Vmax were found at optimal conditions
# # via a paper, and thus no need to extrapolate values.
# activationEnergy = amylaseHelper.GetActivationEnergyFromClausiusClapeyron(
#                                         tempOne=dfToArr[5][5],
#                                         kmOne=dfToArr[5][0],
#                                         tempTwo=dfToArr[9][5],
#                                         kmTwo=dfToArr[9][0])
# reactorTemp = 70+273.15
#
# kmAtReactorConditions = amylaseHelper.GetDesiredKmFromActivationEnergy(
#     kmKnown=dfToArr[5][0], tempKnown=dfToArr[5][5],
#     tempDesired = reactorTemp, activationEnergy=activationEnergy)

'''
MASS BALANCE
'''

initialMassStarch = 5 # g Starch/L
initSubstrateConc = initialMassStarch/359.33 # mol Starch/L

initEnzymeConcentration = np.linspace(1e-7, 9e-7, 5) #M
kcat = constDF['Kcat (1/min)'][9] #1/min
bindingAffinity = constDF['Km (M)'][9] #M
temperature = constDF['Temperature (C)'][9]

residenceTimeArr = np.linspace(0, 480, 500) #min
prodFormedList = []
substrateOutList = []
conversionList = []

masterDataAmylase = pd.DataFrame()
masterDataAmylase['Residence Time'] = residenceTimeArr

reactants = ["Starch"]
products = ["Glucose", "Maltose"]
masterDataAmylase = amylaseHelper.GetConversionVsResidenceTimeWithCstr(initSubstrateConc,
                                                                       initEnzymeConcentration,
                                                                       kcat,
                                                                       bindingAffinity,
                                                                       reactants,
                                                                       products)

titles = masterDataAmylase.columns

fig = plt.figure()

for i in range(4, len(titles), 4):
    plt.plot(masterDataAmylase['Residence Time'], masterDataAmylase[titles[i]], label=titles[i])


plt.xlabel("Residence Time [min]")
plt.ylabel("Conversion [M]")
plt.title("Conversion vs Residence Time - Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, 
                                                                bindingAffinity)
          ,fontweight='bold')
plt.legend(loc="best")

masterDataAmylase.to_csv(Path("Output/Amylase/CSVs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, bindingAffinity))
plt.savefig(Path("Output/Amylase/Graphs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, bindingAffinity))
plt.show()


