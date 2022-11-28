# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:37:44 2022

@author: CRING
"""

from sympy import *

#Dose = dose
#Time = Time
#T_end = end Time
#Assumption: Dose and Time are nonnegative
Dose, Time = symbols('Dose, Time', nonnegative = True)
#Assumption: T_end is positive
T_end = symbols('T_end', positive = True)

# One compartment model analytic AUC
## Parameters:
    ## kelim = Elimination rate
    ## Vdist = Volume of distribution
    ## Assumption: Both of these are positive
kelim, Vdist= symbols('kelim, Vdist', positive = True)
## IV dose
### Cp(Time)
comp1_iv = (Dose/Vdist)*exp(-kelim*Time)
### Integrate from 0 to T_end
integrate(comp1_iv, (Time, 0, T_end))
### Result: Dose/(Vdist*kelim) - Dose*exp(-T_end*kelim)/(Vdist*kelim)

##Oral dose
kgutabs = symbols('kgutabs', positive = True)
Fgutabs = symbols('Fgutabs', nonnegative = True)
comp1_po = (Dose*Fgutabs*kgutabs)/(Vdist*(kgutabs - kelim)) * (exp(-kelim*Time)- exp(-kgutabs*Time))
integrate(comp1_po, (Time, 0, T_end))
### Result: -Dose*Fgutabs*kgutabs*(1/kgutabs - 1/kelim)/(Vdist*(-kelim + kgutabs)) + Dose*Fgutabs*kgutabs*(exp(-T_end*kgutabs)/kgutabs - exp(-T_end*kelim)/kelim)/(Vdist*(-kelim + kgutabs))

### Case when kelim = kgutabs
comp1_po_alt = (Dose*Fgutabs*kelim)/Vdist *Time *exp(-kelim*Time)
integrate(comp1_po_alt, (Time, 0, T_end))
### Result: Dose*Fgutabs/(Vdist*kelim) + (-Dose*Fgutabs*T_end*kelim - Dose*Fgutabs)*exp(-T_end*kelim)/(Vdist*kelim)

## Two compartment model analytic AUC
V1 = symbols('V1', positive = True)
alpha, beta = symbols('alpha, beta', positive = True)
k21, k12 = symbols('k21, k12', positive = True)
A, B = symbols('A, B', positive = True)

### IV data
A_exp = Dose*(alpha - k21)/(V1*(alpha - beta))
B_exp = Dose*(k21 - beta)/(V1*(alpha - beta))
comp2_iv = A*exp(-alpha * Time) + B * exp(-beta*Time)
integrate(comp2_iv, (Time, 0, T_end))
### A/alpha - A*exp(-T_end*alpha)/alpha + B/beta - B*exp(-T_end*beta)/beta

### Oral data
A_exp_po = (kgutabs * Fgutabs/V1 * Dose * (alpha - k21)) / ( (kgutabs - alpha) * (alpha - beta))
B_exp_po = (ka * Fgutabs/V1 * Dose * (k21 - beta)) / ( (kgutabs - beta) * (alpha - beta))

comp2_po =  A * exp(-alpha * Time) + B * exp(-beta * Time) + -(A + B) * exp(-kgutabs * Time)
integrate(comp2_po, (Time, 0, T_end))
### Result: A/alpha - A*exp(-T_end*alpha)/alpha + B/beta - B*exp(-T_end*beta)/beta + (-A - B)/kgutabs - (-A - B)*exp(-T_end*kgutabs)/kgutabs


