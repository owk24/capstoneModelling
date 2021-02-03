# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

invertaseModelHelper = ModelHelper()
invertaseModelHelper.CreateFoldersForOutput("Catalase")
constDF = pd.read_csv(Path("Constants") / "CatalaseConstants.csv")

# initialEnzyme = 9E-09
Kcat = constDF["Kcat (1/min)"][3]
Km = constDF["Km (M)"][3]
# Vm = Kcat * initialEnzyme
temperature = constDF["Temperature (C)"][3]

initialG = 0
initialF = 0
initialS = 10/342.3 #100g/L, 342.3 g/mol molar mass
initialConditions = [initialG, initialF, initialS]


initialEnzymeConcentrations = np.linspace(5e-9, 5e-8, 5)
masterCatalaseOutput = invertaseModelHelper.\
    GetConversionVsResidenceTimeWithCstrForCatalase(initialS, initialEnzymeConcentrations, Kcat, Km)

titles = masterCatalaseOutput.columns
fig = plt.figure()

for i in range(4, len(titles), 5):
    plt.plot(masterCatalaseOutput['Residence Time'], masterCatalaseOutput[titles[i]], label=titles[i])

plt.xlabel("Residence Time [min]")
plt.ylabel("Conversion [M]")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Invertase - Conversion (1:1) vs Residence Time", fontweight='bold')
plt.legend(loc="best")

masterCatalaseOutput.to_csv(Path("Output/Catalase/CSVs") / "Catalase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
plt.savefig(Path("Output/Catalase/Graphs") / "Catalase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()

# t = np.linspace(0, 480, 100)
# invertaseModelOutput = invertaseModelHelper.GetRatesVersusTimeFromModel(t, km, vm, initialConditions)
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
