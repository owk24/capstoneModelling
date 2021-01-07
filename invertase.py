# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 21:32:35 2021

@author: Lakshay Verma and Owen Kim

"""

'''
PACKAGES
Import all relevant packages for analytics
'''

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Helper import modelHelper

'''
Alternative "inhibition" model. To be investigated......
'''

# Ki Km has units of g/L and beta has no units
KiFruc = 47.99
KmFruc = 19.84
betaFruc = 0.421

KiGluc = 10.73
KmGluc = 14.19
betaGluc = 0.64

# Defines the non-competitive inhibitor model that is described in
# https://www.scielo.br/scielo.php?script=sci_arttext&pid=S0104-66321999000200006
def partialNonCompModel(Vm, beta, Ki, Km, S, I):
    return (Vm*S*(1+beta*I)/Ki)/((S+Km)*(1+I/Ki))

# Input conditions for the reaction (g/L)
initialSucroseConc = 10
initialEnzymeConc = 0.1
initialInhibitorConc = 0

constDF = pd.read_csv(Path("Constants") / "InvertaseConstants.csv")
dfToArr = constDF.to_numpy()

invertaseHelper = modelHelper.modelHelper(kmOne=dfToArr[0][0], tempOne=dfToArr[0][5], 
                                          kmTwo=dfToArr[1][0], tempTwo=dfToArr[1][5])

activationEnergy = invertaseHelper.GetActivationEnergyFromClausiusClapeyron()
print(invertaseHelper.GetDesiredKmFromActivationEnergy(invertaseHelper.kmOne, invertaseHelper.tempOne, 36, activationEnergy))


