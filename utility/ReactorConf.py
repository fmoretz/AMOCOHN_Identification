# Reactor Conditions
from PhysConstants import *
import numpy as np

D   = 0.05    # [1/d] Dilution rate
T    = 35      # [°C]
Pt   = 1       # [atm]
alfa = 1

C_d = [0.1, 0.1]

sol  = np.exp(-159.854 + 8741.68/(T+273.15) + 21.6694*np.log(T+273.15) - 0.00110261*(T+273.15))   # [-] mole fraction of dissolved CO2 in water at [T]°C
KH = 1/sol                                                        # [atm] -        Henry's constant CO2 at [T]°C - Partial Pressure Relation
KH = 1/(KH/55342)                                                 # [mmol/L/atm]   Henry's constant CO2 at [T]°C - Concentrations Relation
