# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Helper.ModelHelper import ModelHelper

constDF = pd.read_csv(Path("Constants") / "AmylaseConstants.csv")   
amylaseHelper = ModelHelper()
amylaseHelper.CreateFoldersForOutput("Amylase")

initalMassSucrose = 5 # g/L
initSubstrateConc = initalMassSucrose/359.33 #g/L -> mol/L

initEnzymeConcentration = np.linspace(2e-7, 9e-7, 5) #M
kcat = constDF['Kcat (1/min)'][9]
bindingAffinity = constDF['Km (M)'][9]
temperature = constDF['Temperature (C)'][9]

residenceTimeArr = np.linspace(0, 480, 500) #min
prodFormedList = []
substrateOutList = []
conversionList = []

reactants = ["Starch"]
products = ["Glucose", "Maltose"]
amylaseModelBatchDF = amylaseHelper.GetDataFromBatchModel(bindingAffinity, kcat, initEnzymeConcentration,
                                                          initSubstrateConc, residenceTimeArr, reactants, products)
titles = amylaseModelBatchDF.columns
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(amylaseModelBatchDF["Time [min]"], amylaseModelBatchDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, bindingAffinity), fontweight='bold')
plt.suptitle("Amylase - Conversion vs Time", fontweight='bold')
plt.legend(loc="best")

amylaseModelBatchDF.to_csv(Path("Output/Amylase/CSVs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, bindingAffinity))
plt.savefig(Path("Output/Amylase/Graphs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, bindingAffinity))
plt.show()

