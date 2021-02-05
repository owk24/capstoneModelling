# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

invertaseModelHelper = ModelHelper()
invertaseModelHelper.CreateFoldersForOutput("Invertase")
constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")

listOfReactants = ["Sucrose"]
listOfProducts = ["Glucose", "Fructose"]
Kcat = constDF["Kcat (1/min)"][8]
Km = constDF["Km (M)"][8]
temperature = constDF["Temperature (C)"][8]

initialG = 0 #M
initialF = 0 #M
initialS = 10/342.3 #g/L -> mol/g, 342.3 g/mol molar mass
initialConditions = [initialG, initialF, initialS]

initialEnzymeConcentrations = np.linspace(5e-9, 5e-8, 5)
residenceTimes = np.linspace(0, 480, 480)
invertaseModelBatchDF = invertaseModelHelper.GetDataFromBatchModel(Km, Kcat, initialEnzymeConcentrations, initialS,
                                                                   residenceTimes, ["Sucrose"], ["Glucose", "Fructose"])

invertaseModelBatchDF.to_csv(Path("Output/Invertase/CSVs") / "Invertase_T={0:.0f}C_Km={1:.2e}.csv"
                             .format(temperature, Km))

titles = invertaseModelBatchDF.columns
fig = plt.figure()
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(invertaseModelBatchDF["Time [min]"], invertaseModelBatchDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Invertase - Conversion vs Time", fontweight='bold')
plt.legend(loc="best")

plt.savefig(Path("Output/Invertase/Graphs") / "Invertase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()
