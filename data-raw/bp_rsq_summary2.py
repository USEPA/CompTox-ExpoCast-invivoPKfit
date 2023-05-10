# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:24:59 2023

@author: CRING
"""

from sympy import *
from sympy.abc import i, j
yi = IndexedBase('yi')
mui = IndexedBase('mui')
Ni = IndexedBase('Ni')
b0, b1 = symbols('b0 b1', real = True)

# Numerator
summand = (
            (yi[j] - mui)**2 -
        (b0 + b1 * mui)
        )**2
summand = expand(summand)
summand = collect(summand, yi[j])
mysum = Sum(
        summand,
        (j, 1, Ni)
        )

mysum1 = expand(mysum).doit()

#apply simplify() to each term separately
#get args (terms) as a tuple
myargs = mysum1.args
#convert to a list
myargs2 = list(myargs)
#simplify each element
for i in range(len(myargs2)):
    myargs2[i] = simplify(myargs2[i])
#convert back to tuple
myargs3 = tuple(myargs2)
#reconstruct the expression
mysum2 = Add(*myargs3)

#now substitute
si = symbols('si', positive = True, real = True)
ybari = symbols('ybari', real = True)

mysum3 = mysum2.subs(Sum(yi[j]**2, (j, 1, Ni)), 
                     (Ni-1)*si**2 + Ni * ybari**2)
mysum4 = mysum3.subs(Sum(yi[j], (j, 1, Ni)),
                     Ni * ybari)

#we still have yi[j]**3 and yi[j]**4 sums here
#anything we can do to simplify these?
#if we assume that yi[j] are normally distributed around ybar[i],
#then the third central moment is zero, 
#and the fourth central moment is 3*si**4

N = symbols('N', integer = True, positive = True)
#write the third central moment and equate it to zero
m3 = expand(Sum(expand((yi[j] - ybari)**3), (j, 1, Ni)).doit()/N)
m3args = list(m3.args)
for i in range(len(m3args)):
    m3args[i] = simplify(m3args[i])
#convert back to tuple
m3args = tuple(m3args)
#reconstruct the expression
m3 = Add(*m3args)
#substitute
m3 = m3.subs(Sum(yi[j]**2, (j, 1, Ni)), 
                     (Ni-1)*si**2 + Ni * ybari**2)
m3 = m3.subs(Sum(yi[j], (j, 1, Ni)),
                     Ni * ybari)
#equate it to zero and solve for Sum(yi[j]**3)
m3_exp = simplify(solveset(expand(m3), Sum(yi[j]**3, (j, 1, Ni))))
m3_subexp = m3_exp.args[0] #sub this for Sum(yi[j]**3)

#similarly
m4 = expand(Sum(expand((yi[j] - ybari)**4), (j, 1, Ni)).doit()/N)
m4args = list(m4.args)
for i in range(len(m4args)):
    m4args[i] = simplify(m4args[i])
#convert back to tuple
m4args = tuple(m4args)
#reconstruct the expression
m4 = Add(*m4args)
#substitute
m4 = m4.subs(Sum(yi[j]**3, (j, 1, Ni)),
             m3_subexp)
m4 = m4.subs(Sum(yi[j]**2, (j, 1, Ni)), 
                     (Ni-1)*si**2 + Ni * ybari**2)
m4 = m4.subs(Sum(yi[j], (j, 1, Ni)),
                     Ni * ybari) 

#now equate to 3*si**4 and solve
m4_subexp = simplify(solveset(m4 - 3*si**4, Sum(yi[j]**4, (j, 1, Ni)))).args[0]

#now substitute back into mysum4
mysum5 = mysum4.subs(Sum(yi[j]**3, (j, 1, Ni)),
             m3_subexp)
mysum6 = mysum5.subs(Sum(yi[j]**4, (j, 1, Ni)), m4_subexp) 

mysum7 = collect(collect(collect(collect(expand(mysum6), Ni), 2*b0), 2*b1*mui), ybari)
mysum8 = mysum7.subs(-ybari**2 + 2*ybari*mui - mui**2,
                     -(ybari - mui)**2)
mysum9 = collect(mysum8, 2*si**2)
mysum10 = mysum9.subs(- 3*ybari**2 + 6*ybari*mui - 3*mui**2,
            -3*(ybari - mui)**2)

# Denominator

## Grand average squared residual
#this I can do by hand: it comes to 
grand_avg_eps2_summand = (Ni-1) * si**2 +\
 Ni * ybari**2 -\
 2 * mui * Ni * ybari +\
 Ni*mui**2
 #now make this properly into a sum
si = IndexedBase('si', positive = True, real = True)
ybari = IndexedBase('ybari', real = True)
G = symbols('G', integer = True, positive = True)
grand_avg_eps2 = (1/N) * \
Sum((Ni[i]-1) * si[i]**2 +\
 Ni[i] * ybari[i]**2 -\
 2 * mui[i] * Ni[i] * ybari[i] +\
 Ni[i]*mui[i]**2, (i, 1, G))
    
eps2bar = symbols('eps2bar', real = True, positive = True)
    
denomsummand =  expand(((yi[j] - mui)**2 - eps2bar)**2)
denom = expand(Sum(denomsummand, (j, 1, Ni))).doit() 

#simp0lify term by term
denomargs = list(denom.args)
for i in range(len(denomargs)):
    denomargs[i] = simplify(denomargs[i])
#convert back to tuple
denomargs = tuple(denomargs)
#reconstruct the expression
denom = Add(*denomargs)

#substitute for yi[j]**2, **3, **4 terms as above
denom2 = denom.subs(Sum(yi[j]**4, (j, 1, Ni)), m4_subexp) 
 
denom3 = denom2.subs(Sum(yi[j]**3, (j, 1, Ni)),
             m3_subexp)

denom4 = denom3.subs(Sum(yi[j]**2, (j, 1, Ni)), 
                     (Ni-1)*si**2 + Ni * ybari**2)
denom5 = denom4.subs(Sum(yi[j], (j, 1, Ni)),
                     Ni * ybari) 


