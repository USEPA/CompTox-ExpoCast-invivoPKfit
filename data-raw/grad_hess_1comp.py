# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:39:07 2022

@author: CRING
"""

from sympy import *
import pandas as pd

Dose, Time, logValue, logLOQ = symbols('Dose, Time, logValue, logLOQ', positive = True)
sigma_ref = symbols('sigma_ref', positive = True)
kelim, Vdist = symbols('kelim, Vdist', positive = True)
Fgutabs, kgutabs = symbols('Fgutabs, kgutabs', positive = True)

#1-compartment model gradient & hessian

mu = symbols('mu', positive = True)
z = symbols ('z')
phi, Phi = symbols('phi Phi')
pred = symbols('pred')
x = symbols('x')
 
def grad_log_like(route, wrt, detect):
        if route=='iv':
            #IV 1-compartment model
            pred_exp = Dose*exp(-kelim * Time)/Vdist
            mu_exp = log(pred_exp)
        else:
            #Oral 1-compartment model
            po_const = (Fgutabs * Dose * kgutabs)/((Vdist*(kgutabs - kelim)))
            pred_exp = po_const * (exp(-kelim * Time) - exp(-kgutabs* Time))
            mu_exp = log(pred_exp)
        
        if detect: #PDF
            f = log(1/(sigma_ref*sqrt(2*S.Pi))*exp(-((x - mu)/sigma_ref)**2/2))
        else: #CDF
            f = log((1 + erf((x-mu)/(sigma_ref*sqrt(2)) ))/2)

        if(wrt == sigma_ref):
             df_dvar = diff(f, sigma_ref)
        else:
             #derivative wrt mu
             df_dmu = diff(f, mu)
             #derivative of mu wrt var
             dmu_dvar = diff(mu_exp, wrt)
             #substitute expression for mu
             dmu_dvar = dmu_dvar.subs(mu_exp, mu)
             #substitute in the expression for predicted value
             #dmu_dvar = dmu_dvar.subs(pred_exp, pred)
             #simplify
             #dmu_dvar = simplify(dmu_dvar)
             #chain rule
             df_dvar = df_dmu * dmu_dvar
        
        #simplify
        #subsitute in z
        df_dvar = df_dvar.subs((x - mu)/sigma_ref, z)
        df_dvar = df_dvar.subs((mu - x), -z*sigma_ref)
        #substitute in standard normal PDF, phi,
        #and stdandard normal CDF, Phi
        df_dvar = df_dvar.subs(1/(sigma_ref*sqrt(2*S.Pi))*exp(-z**2/2), phi)
        df_dvar = df_dvar.subs((1 + erf(z/sqrt(2) ))/2, Phi)
        df_dvar = df_dvar.subs(erf(sqrt(2)*z/2) + 1, 2*Phi)
        df_dvar = df_dvar.subs(exp(-z**2/2), sigma_ref*sqrt(2)*sqrt(S.Pi)*phi)
        
        df_dvar = simplify(df_dvar)
        df_dvar = df_dvar.subs((mu - x), -z*sigma_ref)
        df_dvar = simplify(df_dvar)
        return(df_dvar)
    

       
for this_detect in True, False:
    for this_route in 'po', 'iv':
        for this_par in kelim, Vdist, Fgutabs, kgutabs:
            print('For detect = ', this_detect, 'route = ', this_route, 'Gradient of log-likelihood wrt', this_par, ':\n')
            print(grad_log_like(this_route, this_par, this_detect))
            print('\n')
            

print('#### SIGMAS ####\n')
        
for this_detect in True, False:
    for this_route in 'po', 'iv':
            print('For detect = ', this_detect, 'route = ', this_route, 'Gradient of log-likelihood wrt', sigma_ref, ':\n')
            print(grad_log_like(this_route, sigma_ref, this_detect))
            print('\n')        


