# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 11:21:07 2023

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

#deriv wrt time
diff(comp1_iv, Time)

##Oral dose
kgutabs = symbols('kgutabs', positive = True)
Fgutabs = symbols('Fgutabs', nonnegative = True)
comp1_po = (Dose*Fgutabs*kgutabs)/(Vdist*(kgutabs - kelim)) * (exp(-kelim*Time)- exp(-kgutabs*Time))

diff(comp1_po, Time)