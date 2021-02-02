# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

invertaseModelHelper = ModelHelper()
invertaseModelHelper.CreateFoldersForOutput("Invertase")
constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")

initialEnzyme = 9E-09
kcat = constDF["Kcat (1/min)"][8]
km = constDF["Km (M)"][8]
vm = kcat * initialEnzyme
temperature = constDF["Temperature (C)"][8]

initialG = 0
initialF = 0
initialS = 10/342.3 #100g/L, 342.3 g/mol molar mass
initialConditions = [initialG, initialF, initialS]

t = np.linspace(0, 480, 100)
invertaseModelOutput = invertaseModelHelper.GetRatesVersusTimeFromModel(t, km, vm, initialConditions)

exportTitle = "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}".format(temperature, km, kcat)
plt.plot(t, invertaseModelOutput["Glucose [M]"], label="Glucose")
plt.plot(t, invertaseModelOutput["Fructose [M]"], label="Fructose")
plt.plot(t, invertaseModelOutput["Sucrose [M]"], label="Sucrose")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [M]")
plt.title("Temp={0:.2f}C, Km={1:.2e}M, Vm={2:.2e}M/min".format(temperature, km, kcat))
plt.legend(loc="best")
plt.savefig(Path("Output/Invertase/Graphs") /
            "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}.png".format(temperature, km, kcat))
plt.show()
invertaseModelOutput.to_csv(Path("Output/Invertase/CSVs") /
            "Invertase_T={0:.0f}C_Km={1:.2e}_Vm={2:.2e}.csv".format(temperature, km, kcat))
