# Baseline-model
Finds the best baseline parameters for each light curve in CONAN

--------------------
base_model.py
--------------------
Requirements: 
1) install itertools
2) In order to be able to compare the BIC between different fits, the errors must be the same. The CF option (input_v26_BAT.dat) adapts the errors, so make sure it is set to "none". 

Goal: To find the best baseline model for each light curve using the Bayes factor. The code approximates the bayes factor by using different BICs (Bayesian information criterion).
The best model should have the lowest possible BIC and the lowest possible free parameters (npar).

- First the code defines the polynomial orders for each parameter. For time it ranges from 0-4 and for the rest (air mass, xy shifts, FWHM and sky values) it ranges from 0-2.
- Calculates the BIC for all the different baseline parameter combinations (405 combinations) and finds the minimum BIC.
- Calculates the bayes factor using: Bf = exp(ΔBIC/2) and finds the models with a very strong evidence (Bf>150)
- Finally finds the model with the lowest BIC and minimum number of free parameters (the simplest possible model).


* After finding the best baseline parameters, set the CF option to rchisq and update the best baseline parameters in the input_v26_BAT.dat file. *
