# Reactor Conditions
from utility.PhysConstants import *
import numpy as np

T    = 35      # [°C]
Pt   = 1       # [atm]
alpha = 0.7

C_d = [0.1, 0.1]

sol  = np.exp(-159.854 + 8741.68/(T+273.15) + 21.6694*np.log(T+273.15) - 0.00110261*(T+273.15))   # [-] mole fraction of dissolved CO2 in water at [T]°C
KH = 1/sol                                                        # [atm] -        Henry's constant CO2 at [T]°C - Partial Pressure Relation
KH = float(1/(KH/55342))                                                 # [mmol/L/atm]   Henry's constant CO2 at [T]°C - Concentrations Relation
