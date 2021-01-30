# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from Helper.ModelHelper import ModelHelper

invertaseModelHelper = ModelHelper()
constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")

initialEnzyme = 9E-09
kcat = constDF["Kcat (1/min)"][8]
km = constDF["Km (M)"][8]
temperature = constDF["Temperature (C)"][8]

initialG = 0
initialF = 0
initialS = 10/342.3 #100g/L, 342.3 g/mol molar mass
initialConditions = [initialG, initialF, initialS]

t = np.linspace(0, 480, 100)
test = invertaseModelHelper.GetRatesVersusTimeFromModel(t, km, kcat, initialConditions, initialEnzyme)

plt.plot(t, test["Glucose [M/min]"], label="Glucose")
plt.plot(t, test["Fructose [M/min]"], label="Fructose")
plt.plot(t, test["Sucrose [M/min]"], label="Sucrose")
plt.xlabel("Time (min)")
plt.ylabel("Rate (M/min)")
plt.title("Temp={0:.2f}C, Km={1:.2e}M, Vm={2:.2e}M/min"
          .format(temperature, km, kcat))
plt.legend(loc="best")
plt.show()
# concentrationDf.to_csv("./Outputs/Invertase_T={0:.0f}C_E0={1:.2e}M.csv"
#                        .format(constDF["Temperature (C)"][i], 9E-09))
