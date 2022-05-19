import numpy as np

from Identification import*
from Influent import*

mu1_SS = alfa*D + kd[0]
mu2_SS = alfa*D + kd[1]

S1_SS  = Ks[0]*mu1_SS/(mu_max[0] - mu1_SS)

a_coeff = 1/KI2
b_coeff = 1-mu_max[1]/mu2_SS
delta   = b_coeff**2 - 4*Ks[1]*a_coeff

S2_SS   = (-b_coeff - delta**0.5)/2/a_coeff

XT_SS   = XT_in*D/(D+k[6])

X1_SS   = (D*(S1_in - S1_SS) + k[6]*XT_SS)/(k[0]*mu1_SS)
X2_SS   = 1/(k[2]*mu2_SS)*(D*(S2_in - S2_SS) + k[1]/k[0]*(D*(S1_in - S1_SS) + k[6]*XT_SS));

q_M_SS  = k[5]*mu2_SS*X2_SS

gamma   = (k[0]*N_S1 - N_bac)*mu1_SS*X1_SS - N_bac*mu2_SS*X2_SS + kd[0]*N_bac*X1_SS + kd[1]*N_bac*X2_SS
Z_SS    = Z_in + gamma/D

epsi    = D*(C_in + S2_SS - Z_SS) + k[3]*mu1_SS*X1_SS + k[4]*mu2_SS*X2_SS;
om      = Pt*kLa*D*KH + kLa*epsi + q_M_SS*(kLa + D);
Pc_SS   = (om - (om**2-4*kLa**2*epsi*Pt*D*KH)**0.5)/(2*kLa*D*KH);
CO2_SS  = (epsi + kLa*KH*Pc_SS)/(kLa + D);
C_SS    = CO2_SS + Z_SS - S2_SS;
q_C_SS  = kLa*(CO2_SS - KH*Pc_SS);

SSTATE  = [XT_SS, X1_SS, X2_SS, Z_SS, S1_SS, S2_SS, C_SS, CO2_SS, Pc_SS, q_M_SS, q_C_SS];

print('\n STEADY STATE VALUES \n')
print('S.S. of XT  ', float(XT_SS))
print('S.S. of X1  ', float(X1_SS))
print('S.S. of X2  ', float(X2_SS))
print('S.S. of Z   ', float(Z_SS))
print('S.S. of S1  ', float(S1_SS))
print('S.S. of S2  ', float(S2_SS))
print('S.S. of C   ', float(C_SS))
print('S.S. of CO2 ', float(CO2_SS))
print('S.S. of Pc  ', float(Pc_SS))
print('S.S. of q_M ', float(q_M_SS))
print('S.S. of q_C ', float(q_C_SS))
