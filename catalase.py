# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

catalaseModelHelper = ModelHelper()
catalaseModelHelper.CreateFoldersForOutput("Catalase")
constDF = pd.read_csv(Path("Constants") / "CatalaseConstants.csv")

Kcat = constDF["Kcat (1/min)"][3]
Km = constDF["Km (M)"][3]
temperature = constDF["Temperature (C)"][3]

initialG = 0
initialF = 0
initialS = 0.1/34.02 #g/L -<mol/L
initialConditions = [initialG, initialF, initialS]
residenceTimes = np.linspace(0, 480, 500) #min
reactants = ["Hydrogen Peroxide"]
products = ["Water", "Oxygen"]

initialEnzymeConcentrations = np.linspace(9e-11, 5e-10, 5)
catalaseBatchDF = catalaseModelHelper.GetDataFromBatchModel(Km, Kcat, initialEnzymeConcentrations, initialS,
                                                            residenceTimes, reactants, products)

titles = catalaseBatchDF.columns
fig = plt.figure()
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(catalaseBatchDF["Time [min]"], catalaseBatchDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Catalase - Conversion (1:1) vs Time", fontweight='bold')
plt.legend(loc="best")

catalaseBatchDF.to_csv(Path("Output/Catalase/CSVs") / "Catalase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
plt.savefig(Path("Output/Catalase/Graphs") / "Catalase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()

# t = np.linspace(0, 480, 100)
# invertaseModelOutput = catalaseModelHelper.GetRatesVersusTimeFromModel(t, km, vm, initialConditions)
# exportTitle = "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}".format(temperature, Km, Kcat)
# plt.plot(t, invertaseModelOutput["Glucose [M]"], label="Glucose")
# plt.plot(t, invertaseModelOutput["Fructose [M]"], label="Fructose")
# plt.plot(t, invertaseModelOutput["Sucrose [M]"], label="Sucrose")
# plt.xlabel("Time [min]")
# plt.ylabel("Concentration [M]")
# plt.title("Temp={0:.2f}C, Km={1:.2e}M, Vm={2:.2e}M/min".format(temperature, Km, Kcat))
# plt.legend(loc="best")
# plt.savefig(Path("Output/Invertase/Graphs") /
#             "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}.png".format(temperature, Km, Kcat))
# invertaseModelOutput.to_csv(Path("Output/Invertase/CSVs") /
#             "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}.csv".format(temperature, Km, Kcat))
# plt.show()
