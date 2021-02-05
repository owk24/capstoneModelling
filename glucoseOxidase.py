# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Helper.ModelHelper import ModelHelper

constDF = pd.read_csv(Path("Constants") / "GlucoseOxidaseConstants.csv")
glucoseOxidaseHelper = ModelHelper()
glucoseOxidaseHelper.CreateFoldersForOutput("Glucose Oxidase")

initSubstrateConc = 5/180.16 # 5g/L -> M
initEnzymeConcentration = np.linspace(7e-9, 1.5e-8, 5) #M
kcat = constDF['Kcat (1/min)'][4]
bindingAffinity = constDF['Km (M)'][4]
temperature = constDF['Temperature (C)'][4]

residenceTimes = np.linspace(0, 480, 500) #min
prodFormedList = []
substrateOutList = []
conversionList = []

reactants = ["Glucose"]
products = ["Gluconic Acid", "Hydrogen Peroxide"]
GOXBatchModelDF = glucoseOxidaseHelper.GetDataFromBatchModel(bindingAffinity, kcat, initEnzymeConcentration,
                                                             initSubstrateConc, residenceTimes, reactants, products)
titles = GOXBatchModelDF.columns
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(GOXBatchModelDF["Residence Time [min]"], GOXBatchModelDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, bindingAffinity), fontweight='bold')
plt.suptitle("Glucose Oxidase - Conversion vs Residence Time", fontweight='bold')
plt.legend(loc="best")

GOXBatchModelDF.to_csv(Path("Output/Glucose Oxidase/CSVs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, bindingAffinity))
plt.savefig(Path("Output/Glucose Oxidase/Graphs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, bindingAffinity))
plt.show()

