# AMOCO_HN with modified identification

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from SS_Algebraic import* # Import the function to calculate SS
from AMOCO_HN_ODE import f_Model_Deviations_AM2HN

from Influent import*
from Identification import*




# SYSTEM
SSTATE = Steady_States_Assessment(y_in, k, kd, alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, N_bac, N_S1)
y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]]


tspan = np.linspace(1,200,10000)
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

print(f'mu1,max: {mu1_max}; Ks1:  {KS1}; Cd1: {C_d[0]}')
print(f'mu2,max: {mu2_max}; Ks2:  {KS2}; KI2: {KI2}; Cd2: {C_d[1]}')
print("Mole fraction of methane in the gas at the end",float(x_M[-1]))
print("Mass fraction of methane in the gas at the end",float(x_M_W[-1]))
plt.close("all")


plt.figure(1)
ax1 = plt.subplot()
ax1.plot(tspan, X1, label="Acidogenics")
ax1.plot(tspan, X2, label="Methanogens")
ax1.set_xlabel('time [d]')
ax1.set_ylabel('Microbial Concentration [g/L]')
ax1.grid(True)
ax1.legend()

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('pH Value [-]', color = color)
ax2.plot(tspan, pH, linestyle='dashed',color =color, label="pH")
ax2.tick_params(axis='y', labelcolor=color)

plt.tight_layout()

plt.figure(2)
sub1 = plt.subplot(2,2,1)
sub1.plot(tspan,q_M,label="CH4")
sub1.plot(tspan,q_C,label="CO2")
sub1.set_ylabel('Gas Flowrate [mmol/L/d]')
sub1.set_xlabel('Time [d]')
sub1.grid(True)
sub1.legend()

sub2 = plt.subplot(2,2,2)
sub2.plot(tspan,Z,label="Alkalinity")
sub2.plot(tspan,C,label="In. Carbon")
sub2.set_ylabel('Inorganics Conc. [mmol/L]')
sub2.set_xlabel('Time [d]')
sub2.grid(True)
sub2.legend()

sub3 = plt.subplot(2,2,3)
sub3.plot(tspan,XT,label="Particulate (XT)")
sub3.plot(tspan,S1,label="COD (S1)")
sub3.set_ylabel('Substrates Conc. [g/L]')
sub3.set_xlabel('Time [d]')
sub3.grid(True)
sub3.legend()

sub4 = plt.subplot(2,2,4)
sub4.plot(tspan,S1,label="COD (S1)")
sub4.plot(tspan,S2,label="VFA (S2) [mmol/L]")
sub4.set_ylabel('Substrates Conc. ')
sub4.set_xlabel('Time [d]')
sub4.grid(True)
sub4.legend()

plt.figure(3)
sub1 = plt.subplot(6,2,1)
sub1.plot(tspan,S1/S1[0])
sub1.set_ylabel('S1*')
sub1.tick_params(labelbottom=False)
sub1.set_xlim(tspan[0],tspan[-1])
sub1.grid(True)

sub2 = plt.subplot(6,2,2)
sub2.plot(tspan,S2/S2[0])
sub2.set_ylabel('S2*')
sub2.tick_params(labelbottom=False)
sub2.set_xlim(tspan[0],tspan[-1])
sub2.grid(True)

sub3 = plt.subplot(6,2,3)
sub3.plot(tspan,X1/X1[0])
sub3.set_ylabel('X1*')
sub3.tick_params(labelbottom=False)
sub3.set_xlim(tspan[0],tspan[-1])
sub3.grid(True)

sub4 = plt.subplot(6,2,4)
sub4.plot(tspan,X2/X2[0])
sub4.set_ylabel('X2*')
sub4.tick_params(labelbottom=False)
sub4.set_xlim(tspan[0],tspan[-1])
sub4.grid(True)

sub5 = plt.subplot(6,2,5)
sub5.plot(tspan,C/C[0])
sub5.set_ylabel('C*')
sub5.tick_params(labelbottom=False)
sub5.set_xlim(tspan[0],tspan[-1])
sub5.grid(True)

sub6 = plt.subplot(6,2,6)
sub6.plot(tspan,Z/Z[0])
sub6.set_ylabel('Z*')
sub6.tick_params(labelbottom=False)
sub6.set_xlim(tspan[0],tspan[-1])
sub6.grid(True)

sub7 = plt.subplot(6,2,7)
sub7.plot(tspan,CO2/CO2[0])
sub7.set_ylabel('CO2*')
sub7.tick_params(labelbottom=False)
sub7.set_xlim(tspan[0],tspan[-1])
sub7.grid(True)

sub8 = plt.subplot(6,2,8)
sub8.plot(tspan,B/B[0])
sub8.set_ylabel('B*')
sub8.tick_params(labelbottom=False)
sub8.set_xlim(tspan[0],tspan[-1])
sub8.grid(True)

sub9 = plt.subplot(6,2,9)
sub9.plot(tspan,pH/pH[0])
sub9.set_ylabel('pH*')
sub9.tick_params(labelbottom=False)
sub9.set_xlim(tspan[0],tspan[-1])
sub9.grid(True)

sub10 = plt.subplot(6,2,10)
sub10.plot(tspan,XT/XT[0])
sub10.set_ylabel('XT*')
sub10.tick_params(labelbottom=False)
sub10.set_xlim(tspan[0],tspan[-1])
sub10.grid(True)

sub11 = plt.subplot(6,2,11)
sub11.plot(tspan,q_C/q_C[0])
sub11.set_ylabel('q_C*')
sub11.set_xlabel('Time [d]')
sub11.set_xlim(tspan[0],tspan[-1])
sub11.grid(True)

sub12 = plt.subplot(6,2,12)
sub12.plot(tspan,q_M/q_M[0])
sub12.set_ylabel('q_M*')
sub12.set_xlabel('Time [d]')
sub12.set_xlim(tspan[0],tspan[-1])
sub12.grid(True)

# plt.show()
