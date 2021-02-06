# -*- coding: utf-8 -*-
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Helper.ModelHelper import ModelHelper

amylaseHelper = ModelHelper()
amylaseHelper.CreateFoldersForOutput("Amylase")
inputsDF = pd.read_csv(Path("./") / "AmylaseInputs.csv")
constDF = pd.read_csv(amylaseHelper.ResourcePath(Path("Constants") / "AmylaseConstants.csv"))

Km = float(inputsDF["Km [M]"])
Kcat = float(inputsDF["Kcat [1/s]"])
temperature = int(inputsDF["Temperature [C]"])
initSubstrateConcentration = float(inputsDF["InitalSubstrateConcentration [M]"]) #g/L -> mol/L
initialEnzymeConcentration = np.linspace(inputsDF["InitialEnzymeConcentration [M]"][0],
                                      inputsDF["FinalEnzymeConcentration [M]"][0],
                                      int(inputsDF["AmountOfEnzymeDataPoints"])) #M
residenceTimes = np.linspace(int(inputsDF["InitialTime [min]"]),
                             int(inputsDF["FinalTime [min]"]),
                             int(inputsDF["AmountOfTimePoints"])) #min

reactants = ["Starch"]
products = ["Glucose", "Maltose"]
amylaseModelBatchDF = amylaseHelper.GetDataFromBatchModel(Km, Kcat, initialEnzymeConcentration,
                                                          initSubstrateConcentration, residenceTimes,
                                                          reactants, products)
titles = amylaseModelBatchDF.columns
for i in range(4, len(titles), 4):
    title = titles[i]
    plt.plot(amylaseModelBatchDF["Time [min]"], amylaseModelBatchDF[title], label=title)

plt.xlabel(titles[0])
plt.ylabel("Conversion (%)")
plt.title("Temp={0:.2f}($^\circ$C), Km={1:.2e}M".format(temperature, Km), fontweight='bold')
plt.suptitle("Amylase - Conversion vs Time", fontweight='bold')
plt.legend(loc="best")

amylaseModelBatchDF.to_csv(Path("Output/Amylase/CSVs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.csv".format(temperature, Km))
plt.savefig(Path("Output/Amylase/Graphs") /
            "Amylase_T={0:.0f}C_Km={1:.2e}.png".format(temperature, Km))
plt.show()
