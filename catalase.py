# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper.ModelHelper import ModelHelper

catalaseModelHelper = ModelHelper()
catalaseModelHelper.CreateFoldersForOutput("Catalase")
inputsDF = pd.read_csv(Path("./") / "CatalaseInputs.csv")
constDF = pd.read_csv(catalaseModelHelper.ResourcePath(Path("Constants") / "CatalaseConstants.csv"))

Km = float(inputsDF["Km [M]"])
Kcat = float(inputsDF["Kcat [1/s]"])
temperature = int(inputsDF["Temperature [C]"])
initialSubstrateConcentration = float(inputsDF["InitalSubstrateConcentration [M]"]) #g/L -> mol/L
initialEnzymeConcentrations = np.linspace(inputsDF["InitialEnzymeConcentration [M]"][0],
                                          inputsDF["FinalEnzymeConcentration [M]"][0],
                                          int(inputsDF["AmountOfEnzymeDataPoints"])) #M
residenceTimes = np.linspace(int(inputsDF["InitialTime [min]"]),
                             int(inputsDF["FinalTime [min]"]),
                             int(inputsDF["AmountOfTimePoints"])) #min

reactants = ["Hydrogen Peroxide"]
products = ["Water", "Oxygen"]

catalaseModelBatchDF = catalaseModelHelper.GetDataFromBatchModel(Km, Kcat, initialEnzymeConcentrations,
                                                                 initialSubstrateConcentration, residenceTimes,
                                                                 reactants, products)
titles = catalaseModelBatchDF.columns
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(catalaseModelBatchDF["Time [min]"], catalaseModelBatchDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Catalase - Conversion (1:1) vs Time", fontweight='bold')
plt.legend(loc="best")

catalaseModelBatchDF.to_csv(Path("Output/Catalase/CSVs") /
                       "Catalase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
plt.savefig(Path("Output/Catalase/Graphs") /
            "Catalase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()
