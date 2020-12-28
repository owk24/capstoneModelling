# Project Description
Modelling various enzymatic reactions for the purpose of creating honey without the use of bees. This is a
Capstone Project for the University of Waterloo - Department of Chemical Engineering.

The reactions being modelled are:
* `Sucrose` + `Water` into `Glucose` + `Fructose` with `Invertase`
* `Starch` into `Maltose` + `Fructose` with `Amylase`
* `Glucose` + `Water` + `Oxygen` into `Hydrogen Peroxide` + `Gluconic Acid` with `Glucose Oxidase`
* `Hydrogen Peroxide` into `Water` + `Oxygen` with `Catalase`

### Background Information
Given that this is a design project, Michaelis-Menton kinetics will be used to predict the rate of the reactions.
This is for ease of use and simplicity because there is no access to labs to determine the kinetic properties
of the various enzymes.

# Technical Information
Python version used is `3.8.5`

Packages used are listed in requirements.txt
To avoid conflicts, venv is recommended. Instructions on how to use venv is listed
[here](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/).
**Make sure to install the packages when you have the venv activated and not before**

To install the packages listed in requirements.txt: `pip install -r requirements.txt` while in the same directory as
the txt file.

If Pycharm is being used, [this link](https://www.jetbrains.com/help/pycharm/creating-virtual-environment.html) is
helpful in configuring the environment needed to run the python files with venv.

Members: Owen Kim and Lakshay Verma

