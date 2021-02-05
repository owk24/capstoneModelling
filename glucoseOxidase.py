# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 18:21:07 2021

@author: Lakshay Verma
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

constDF = pd.read_csv(Path("Constants") / "GlucoseOxidaseConstants.csv")   
dfToArr = constDF.to_numpy()

glucoseOxidaseHelper = ModelHelper()
glucoseOxidaseHelper.CreateFoldersForOutput("Glucose Oxidase")

'''
MASS BALANCE
'''


initialMassStarch = 5 # g glucose/L
initSubstrateConc = initialMassStarch/180.16 # mol glucose/L

initEnzymeConcentration = np.linspace(7e-9, 1.5e-8, 5) #M
kcat = constDF['Kcat (1/min)'][4] #1/min
bindingAffinity = constDF['Km (M)'][4] #M
temperature = constDF['Temperature (C)'][4]

residenceTimeArr = np.linspace(0, 480, 500) #min
prodFormedList = []
substrateOutList = []
conversionList = []

masterDataGlucoseOxidase = pd.DataFrame()

reactants = ["Glucose"]
products = ["Gluconic Acid", "Hydrogen Peroxide"]
masterDataGlucoseOxidasecCSTRVOID = glucoseOxidaseHelper.GetConversionVsResidenceTimeWithCstr(initSubstrateConc,
                                                                       initEnzymeConcentration,
                                                                       kcat,
                                                                       bindingAffinity,
                                                                       reactants,
                                                                       products)
masterDataGlucoseOxidase = glucoseOxidaseHelper.GetDataFromBatchModel(bindingAffinity, kcat, initEnzymeConcentration, residenceTimeArr,
                              reactants, products, initSubstrateConc)
masterDataGlucoseOxidase['Residence Time'] = residenceTimeArr
firstColumn = masterDataGlucoseOxidase.pop("Residence Time")
masterDataGlucoseOxidase.insert(0, "Residence Time", firstColumn)

titles = masterDataGlucoseOxidase.columns

fig = plt.figure()

for i in range(4, len(titles), 4):
    plt.plot(masterDataGlucoseOxidase['Residence Time'], masterDataGlucoseOxidase[titles[i]], label=titles[i])


plt.xlabel("Residence Time [min]")
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, 
                                                                bindingAffinity)
          ,fontweight='bold')
plt.suptitle("Glucose Oxidase - Conversion vs Residence Time", fontweight='bold')
plt.legend(loc="best")

masterDataGlucoseOxidase.to_csv(Path("Output/Glucose Oxidase/CSVs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, bindingAffinity))
plt.savefig(Path("Output/Glucose Oxidase/Graphs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, bindingAffinity))


