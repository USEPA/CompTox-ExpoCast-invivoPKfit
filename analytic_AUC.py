# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:37:44 2022

@author: CRING
"""

from sympy import *

#D = dose
#t = time
#T = end time
#Assumption: D and t are nonnegative
D, t = symbols('D, t', nonnegative = True)
#Assumption: T is positive
T = symbols('T', positive = True)

# One compartment model analytic AUC
## Parameters:
    ## ke = Elimination rate
    ## Vd = Volume of distribution
    ## Assumption: Both of these are positive
ke, Vd= symbols('ke, Vd', positive = True)
## IV dose
### Cp(t)
comp1_iv = (D/Vd)*exp(-ke*t)
### Integrate from 0 to T
integrate(comp1_iv, (t, 0, T))
### Result:D/(Vd*ke) - D*exp(-T*ke)/(Vd*ke)


##Oral dose
ka = symbols('ka', positive = True)
Fa = symbols('Fa', nonnegative = True)
comp1_po = (D*Fa*ka)/(Vd*(ka - ke)) * (exp(-ke*t)- exp(-ka*t))
integrate(comp1_po, (t, 0, T))
### Result: -D*Fa*ka*(-1/ke + 1/ka)/(Vd*(ka - ke)) + D*Fa*ka*(-exp(-T*ke)/ke + exp(-T*ka)/ka)/(Vd*(ka - ke))

### Case when kelim = ka
comp1_po_alt = (D*Fa*ke)/Vd *t *exp(-ke*t)
integrate(comp1_po_alt, (t, 0, T))
### Result: D*Fa/(Vd*ke) + (-D*Fa*T*ke - D*Fa)*exp(-T*ke)/(Vd*ke)

## Two compartment model analytic AUC
V1 = symbols('V1', positive = True)
alpha, beta = symbols('alpha, beta', positive = True)
k21, k12 = symbols('k21, k12', positive = True)
A, B = symbols('A, B', positive = True)
A_exp = D*(alpha - k21)/(V1*(alpha - beta))
B_exp = D*(k21 - beta)/(V1*(alpha - beta))
comp2_iv = A*exp(-alpha * t) + B * exp(-beta*t)
integrate(comp2_iv, (t, 0, T))
### A/alpha - A*exp(-T*alpha)/alpha + B/beta - B*exp(-T*beta)/beta

comp2_po = A * exp(-alpha * t) + B * exp(-beta * t) + -(A+B) * exp(-ka * t)
integrate(comp2_po, (t, 0, T))
### Result:A/alpha - A*exp(-T*alpha)/alpha + B/beta - B*exp(-T*beta)/beta + (-A - B)/ka - (-A - B)*exp(-T*ka)/ka


