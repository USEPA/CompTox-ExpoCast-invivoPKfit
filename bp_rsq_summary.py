# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 17:35:31 2023

@author: CRING
"""

from sympy import *
from sympy.abc import i, j

mui = symbols('mui', nonnegative = True, real = True)
dsqbar = symbols('dsqbar', nonnegative = True, real = True)
mubar = symbols('mubar', nonnegative = True, real = True)
Ni = symbols('Ni', positive = True, integer = True, finite = True)
#dsqi = IndexedBase('dsqi')
yi = IndexedBase('yi')


dsqi = expand((yi[j] - mui)*(yi[j] - mui))
summand = expand((dsqi - dsqbar)*(mui - mubar))
mysum = Sum(summand, (j, 1, Ni)).doit()
mysum2 = expand(mysum).doit()

ybari = symbols('ybari', nonnegative = True, real = True)
s2i = symbols('s2i', nonnegative = True, real = True)
mysum3 = mysum2.subs(Sum(-mubar*yi[j]**2, (j, 1, Ni)), mubar * ((Ni-1) * s2i + Ni * ybari**2))
mysum4 = mysum3.subs(Sum(mui*yi[j]**2, (j, 1, Ni)), mui * ((Ni-1) * s2i + Ni * ybari**2))
mysum5 = mysum4.subs(Sum(-2*mui**2*yi[j], (j, 1, Ni)), -2*mui**2*Ni*ybari)
mysum6 = mysum5.subs(Sum(2*mubar*mui*yi[j], (j, 1, Ni)), 2 * mubar * mui * Ni * ybari)
mysum7 = expand(mysum6) #cancellation occurs

N = symbols('N', positive = True, integer = True, finite = True)
dsqbar_exp = expand((1/N) * Sum(dsqi, (j, 1, Ni))).doit()
dsqbar_exp2 = dsqbar_exp.subs(Sum(-2*mui*yi[j], (j, 1, Ni)), -2*mui*Ni*ybari)
dsqbar_exp3 = dsqbar_exp2.subs(Sum(yi[j]**2, (j, 1, Ni)), (Ni-1)*s2i + Ni*ybari**2)
mysum8 = mysum7.subs(dsqbar, dsqbar_exp3)


mubar_exp = (1/N) * Ni * mui
mysum9 = mysum8.subs(mubar, mubar_exp)
mysum10 = expand(mysum9)
mysum11 = collect(mysum10, N)
mysum12 = collect(mysum11, Ni*mui)
mysum13 = collect(mysum12, Ni)
mysum14 = collect(mysum12, mui)