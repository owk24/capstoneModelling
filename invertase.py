# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

invertaseModelHelper = ModelHelper()
invertaseModelHelper.CreateFoldersForOutput("Invertase")
constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")

# initialEnzyme = 9E-09
listOfReactants = ["Sucrose"]
listOfProducts = ["Glucose", "Fructose"]
Kcat = constDF["Kcat (1/min)"][8]
Km = constDF["Km (M)"][8]
# Vm = Kcat * initialEnzyme
temperature = constDF["Temperature (C)"][8]

initialG = 0
initialF = 0
initialS = 10/342.3 #100g/L, 342.3 g/mol molar mass
initialConditions = [initialG, initialF, initialS]


initialEnzymeConcentrations = np.linspace(5e-9, 5e-8, 5)
residenceTimes = np.linspace(0, 480, 480)
masterInvertaseOutput = invertaseModelHelper.GetDataFromBatchModel \
                        (Km, Kcat, initialEnzymeConcentrations, residenceTimes,
                         ["Sucrose"], ["Glucose", "Fructose"])
masterInvertaseOutput.to_csv(Path("Output/Invertase/CSVs") / "Invertase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
# titles = masterInvertaseOutput.columns
# fig = plt.figure()
#
#
# plt.xlabel("Residence Time [min]")
# plt.ylabel("Conversion [M]")
# plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
# plt.suptitle("Invertase - Conversion vs Residence Time", fontweight='bold')
# plt.legend(loc="best")
#
#
# plt.savefig(Path("Output/Invertase/Graphs") / "Invertase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
# plt.show()

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
