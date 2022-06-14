''' Modelling the AMOCO_HN equation in two cases:
    1. Original parameters from the "A generic and systematic procedure to derive a simplified model from
       the anaerobic digestion model No. 1 (ADM1)"; Hassam S. et al., 2015
    2. Modified parameters with the new identification method proposed
    Both cases modeled as a CSTR (alpha = 1)'''

import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve

from AMOCO_HN_ODE import f_Model_Deviations_AM2HN
from Identification import *
from Influent import *
from SS_Algebraic import *  # Import the function to calculate SS

# SYSTEM
tspan = np.linspace(1,160,10000)

''' Case 1: Original parameters'''
# Parameters
k_or = [20, 464, 514, 310, 600, 253, 5.02]
mu_max_or = [0.33, 0.13]
Ks_or = [0.4, 2.93]
KH_or = 16
KI2_or = 207
kLa_or = 24
kd_or = [0.1*mu_max_or[0], 0.1*mu_max_or[1]]
SSTATE = Steady_States_Assessment(y_in, k_or, kd_or, alfa, mu_max_or, Ks_or, KI2_or, KH_or, Pt, kLa_or, D, N_bac, N_S1)
y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]]

YOUT_or = odeint(f_Model_Deviations_AM2HN,y0,tspan,args=(alfa,mu_max_or,Ks_or,KI2_or,KH_or,Pt,kLa_or,D,y_in,k_or,kd_or,N_bac,N_S1))
alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,y_in,k,kd,N_bac,N_S1
XT_or = YOUT_or[:,0]              # [gCOD/L] - Particulate 
X1_or = YOUT_or[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2_or = YOUT_or[:,2]              # [g/L]    - Methanogenic Bacteria
Z_or  = YOUT_or[:,3]              # [mmol/L] - Total Alkalinity
S1_or = YOUT_or[:,4]              # [g/L]    - Organic Soluble Substrate
S2_or = YOUT_or[:,5]              # [mmol/L] - VFA dissolved
C_or  = YOUT_or[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

# Solver Output
mu1_or = np.empty(len(XT_or))
mu2_or = np.empty(len(XT_or))
CO2_or = np.empty(len(XT_or))
B_or   = np.empty(len(XT_or))
phi_or = np.empty(len(XT_or))
p_C_or = np.empty(len(XT_or))
q_C_or = np.empty(len(XT_or))
q_M_or = np.empty(len(XT_or))
pH_or  = np.empty(len(XT_or))


for x in range(len(XT_or)):
    mu1_or[x] = mu_max_or[0]*(S1_or[x]/(S1_or[x]+Ks_or[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2_or[x] = mu_max_or[1]*(S2_or[x]/(S2_or[x]+Ks_or[1]+S2_or[x]**2/KI2_or))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2_or[x] = C_or[x] + S2_or[x] - Z_or[x]                                 # [mmol/L]   - Dissolved CO2
    B_or[x]   = Z_or[x] - S2_or[x]                                        # [mmol/L]   - Alkalinity
    phi_or[x] = CO2_or[x] + KH_or*Pt + k_or[5]/kLa_or*mu2_or[x]*X2_or[x]
    p_C_or[x]  = (phi_or[x] - (phi_or[x]**2- 4*KH_or*Pt*CO2_or[x])**0.5)/(2*KH_or) # [atm]      - CO2 Partial Pressure
    q_C_or[x] = kLa_or*(CO2_or[x] - KH_or*p_C_or[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M_or[x] = k_or[5]*mu2_or[x]*X2_or[x]                                   # [mmol/L/d] - CH4 Outlet Molar Flow
    pH_or[x]  = np.real(-np.log10(Kb*CO2_or[x]/B_or[x]))                  # [-]        - System pH

q_tot_or = q_C_or+q_M_or                                                  # [mmol/L/d] - Outlet global molar flow  
x_M_or   = np.divide(q_M_or,q_tot_or)                                     # [-]        - CH4 Mole fraction
q_C_W_or = q_C_or*44/1000                                              # [g/L/d]    - CO2 Outlet mass flow of
q_M_W_or = q_M_or*16/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W_or = q_C_W_or + q_M_W_or                                          # [g/L/d]    - Outlet global mass flow  
x_M_W_or   = q_M_W_or/q_tot_W_or                                          # [-]        - CH4 Weight Fraction      


print('*** Case 1. Original Parameters***')
print(f"Mole fraction of methane in the gas at the end: {x_M_or[-1]}")
print(f"Mass fraction of methane in the gas at the end: {x_M_W_or[-1]}")

SSTATE = Steady_States_Assessment(y_in, k, kd, alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, N_bac, N_S1)
y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]]

YOUT = odeint(f_Model_Deviations_AM2HN,y0,tspan,args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, y_in, k, kd, N_bac, N_S1))
XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

# Solver Output
mu1 = np.empty(len(XT))
mu2 = np.empty(len(XT))
CO2 = np.empty(len(XT))
B   = np.empty(len(XT))
phi = np.empty(len(XT))
p_C = np.empty(len(XT))
q_C = np.empty(len(XT))
q_M = np.empty(len(XT))
pH  = np.empty(len(XT))



for x in range(len(XT)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2[x] = C[x] + S2[x] - Z[x]                                 # [mmol/L]   - Dissolved CO2
    B[x]   = Z[x] - S2[x]                                        # [mmol/L]   - Alkalinity
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH) # [atm]      - CO2 Partial Pressure
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M[x] = k[5]*mu2[x]*X2[x]                                   # [mmol/L/d] - CH4 Outlet Molar Flow
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))                  # [-]        - System pH

q_tot = q_C+q_M                                                  # [mmol/L/d] - Outlet global molar flow  
x_M   = np.divide(q_M,q_tot)                                     # [-]        - CH4 Mole fraction
q_C_W = q_C*44/1000                                              # [g/L/d]    - CO2 Outlet mass flow of
q_M_W = q_M*16/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W = q_C_W + q_M_W                                          # [g/L/d]    - Outlet global mass flow  
x_M_W   = q_M_W/q_tot_W                                          # [-]        - CH4 Weight Fraction      

print('*** Case 2. New Parameters:***')
print("Mole fraction of methane in the gas at the end",float(x_M[-1]))
print("Mass fraction of methane in the gas at the end",float(x_M_W[-1]))




      