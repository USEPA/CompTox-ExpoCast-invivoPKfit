# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:59:36 2022

@author: CRING
"""
from sympy import *

x, mu, sigma = symbols('x mu sigma')

# Gradient of log of normal PDF

f = log(1/(sigma*sqrt(2*S.Pi))*exp(-((x - mu)/sigma)**2/2))

df_dmu = simplify(diff(f, mu))

df_dsigma = diff(f, sigma)

z = symbols ('z')

df_dmu2 = df_dmu.subs((x - mu)/sigma, z)

df_dsigma2 = df_dsigma.subs((x-mu)/sigma, z)
df_dsigma3 = df_dsigma2.simplify()

#Gradient of log of normal CDF

F = log((1 + erf((x-mu)/(sigma*sqrt(2)) ))/2)

dF_dmu = diff(F, mu)
dF_dmu2 = dF_dmu.subs((x-mu)/sigma, z)

#substitute in standard normal PDF, phi,
#and stdandard normal CDF, Phi
phi, Phi = symbols('phi Phi')
dF_dmu3 = dF_dmu2.subs(1/(sigma*sqrt(2*S.Pi))*exp(-z**2/2), phi)
dF_dmu4 = dF_dmu3.subs((1 + erf(z/sqrt(2) ))/2, Phi)

dF_dsigma = diff(F, sigma)
dF_dsigma2 = dF_dsigma.subs((x-mu)/sigma, z)
dF_dsigma3 = dF_dsigma2.subs(1/(sigma*sqrt(2*S.Pi))*exp(-z**2/2), phi)
dF_dsigma4 = dF_dsigma3.subs((1 + erf(z/sqrt(2) ))/2, Phi)

#Final answers
print('Jacobian\n')
print('Jacobian for detects:\n')
print('df/dmu:', df_dmu2, '\n')
print('df/dsigma:', df_dsigma3, '\n')

print('Jacobian for nondetects:\n')
print('dF/dmu:', dF_dmu4, '\n')
print('dF/dsigma:', dF_dsigma4, '\n')

## Second derivatives??
print('\n\n')
print('Hessian\n')
print('Hessian for detects:\n')

d2f_dmu2 = diff(f, mu, 2)
print('d2f/dmu2:', d2f_dmu2, '\n')

d2f_dmu_dsigma = diff(f, mu, sigma)
print('d2f/dmu_dsigma:', d2f_dmu_dsigma.subs((mu-x), -z*sigma), '\n')


d2f_dsigma_dmu = diff(f, sigma, mu)
print('d2f/dsigma_dmu:', simplify(d2f_dsigma_dmu.subs((mu-x), -z*sigma)), '\n')

d2f_dsigma2 = diff(f, sigma, 2)
print('d2f/dsigma2:', simplify(d2f_dsigma2.subs((mu-x), -z*sigma)), '\n')

print('Hessian for non-detects:\n')
d2F_dmu2 = diff(F, mu, 2)
#substitutions & simplification
d2F_dmu2_s1 = d2F_dmu2.subs((mu-x), -z*sigma)
d2F_dmu2_s2 = d2F_dmu2_s1.subs((-erf(sqrt(2)*z/2) - 1), -2*Phi)
d2F_dmu2_s3 = d2F_dmu2_s2.subs(exp(-z**2/2), sigma*sqrt(2)*sqrt(S.Pi)*phi)
print('d2F/dmu2:', simplify(d2F_dmu2_s3), '\n')

d2F_dmu_dsigma = diff(F, mu, sigma)
d2F_dmu_dsigma_s1 = d2F_dmu_dsigma.subs((mu-x), -z*sigma)
d2F_dmu_dsigma_s2 = d2F_dmu_dsigma_s1.subs((-erf(sqrt(2)*z/2) - 1), -2*Phi)
d2F_dmu_dsigma_s3 = d2F_dmu_dsigma_s2.subs(exp(-z**2/2), sigma*sqrt(2)*sqrt(S.Pi)*phi)
print('d2F/dmu_dsigma:', simplify(d2F_dmu_dsigma_s3), '\n')

d2F_dsigma_dmu = diff(F, sigma, mu)
d2F_dsigma_dmu_s1 = d2F_dsigma_dmu.subs((mu-x), -z*sigma)
d2F_dsigma_dmu_s2 = d2F_dsigma_dmu_s1.subs((-erf(sqrt(2)*z/2) - 1), -2*Phi)
d2F_dsigma_dmu_s3 = d2F_dsigma_dmu_s2.subs(exp(-z**2/2), sigma*sqrt(2)*sqrt(S.Pi)*phi)
print('d2F/dsigma_dmu:', simplify(d2F_dsigma_dmu_s3), '\n')

d2F_dsigma2 = diff(F, sigma, 2)
d2F_dsigma2_s1 = d2F_dsigma2.subs((mu-x), -z*sigma)
d2F_dsigma2_s2 = d2F_dsigma2_s1.subs((-erf(sqrt(2)*z/2) - 1), -2*Phi)
d2F_dsigma2_s3 = d2F_dsigma2_s2.subs(exp(-z**2/2), sigma*sqrt(2)*sqrt(S.Pi)*phi)
print('d2F/dsigma2:', simplify(d2F_dsigma2_s3), '\n')
