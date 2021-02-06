# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Helper.ModelHelper import ModelHelper

glucoseOxidaseHelper = ModelHelper()
glucoseOxidaseHelper.CreateFoldersForOutput("Glucose Oxidase")
inputsDF = pd.read_csv(Path("./") / "GlucoseOxidaseInputs.csv")
constDF = pd.read_csv(glucoseOxidaseHelper.ResourcePath(Path("Constants") / "GlucoseOxidaseConstants.csv"))

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


reactants = ["Glucose"]
products = ["Gluconic Acid", "Hydrogen Peroxide"]

GOXBatchModelDF = glucoseOxidaseHelper.GetDataFromBatchModel(Km, Kcat, initialEnzymeConcentrations,
                                                             initialSubstrateConcentration, residenceTimes,
                                                             reactants, products)
titles = GOXBatchModelDF.columns
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(GOXBatchModelDF["Time [min]"], GOXBatchModelDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Glucose Oxidase - Conversion vs Time", fontweight='bold')
plt.legend(loc="best")

GOXBatchModelDF.to_csv(Path("Output/Glucose Oxidase/CSVs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
plt.savefig(Path("Output/Glucose Oxidase/Graphs") /
            "Gluocse Oxidase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()
