# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from Helper import modelHelper

concentrationOfComponents = []
constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")
kmOne = constDF["Km (M)"][3]
tempOne = constDF["Temperature (K)"][3]
kmTwo = constDF["Km (M)"][4]
tempTwo = constDF["Temperature (K)"][4]

invertaseHelper = modelHelper.modelHelper(kmOne, tempOne, kmTwo, tempTwo)

rangeOfKm = constDF["Km (M)"]
rangeOfVm = constDF["Vm (M/min)"]

def model(inputs, t, index) -> [float, float, float]:
    Vm = constDF["Vm (M/min)"][index]
    Km = constDF["Km (M)"][index]
    r = Vm*inputs[2]/(Km + inputs[2])
    dGdt = r
    dFdt = r
    dSdt = -r
    outputs = [dGdt, dFdt, dSdt]
    return outputs

initialG = 0
initialF = 0
initialS = 10/342.3 #100g/L, 342.3 g/mol molar mass

initialConditions = [initialG, initialF, initialS]

t = np.linspace(0, 60, 100)
i = len(rangeOfKm)-1
concentrationOfComponents = odeint(model, initialConditions, t, args=(i,))

plt.plot(t, concentrationOfComponents[:, 0], label="Glucose")
plt.plot(t, concentrationOfComponents[:, 1], label="Fructose")
plt.plot(t, concentrationOfComponents[:, 2], label="Sucrose")
plt.xlabel("Time (min)")
plt.ylabel("Concentration (M)")
plt.title("Temp={0:.2f} C, Km={1:.2e} M, Vm={2:.2e} M/min"
          .format(constDF["Temperature (C)"][i], constDF["Km (M)"][i], constDF["Vm (M/min)"][i]))
plt.legend(loc="best")
plt.show()

concentrationDf = pd.DataFrame({"Glucose": concentrationOfComponents[:, 0],
                                "Fructose": concentrationOfComponents[:, 1],
                                "Sucrose": concentrationOfComponents[:, 2]})
print(concentrationDf["Glucose"].iloc[-1]*342.3)
