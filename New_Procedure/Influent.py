import numpy as np

from PhysConstants import*
from dataimport import*

# Substrates

N_bac   = 0.08/14*1000*1.55    # Nitrogen content of the biomass [mmolN/gVS]
N_aa    = 0.007*1000*1.55      # Nitrogen content of aminoacids and proteins [mmolN/gVS]
f_aa    = 0.00003              # fraction of aminoacids in the organics substrate S1
f_pr    = 0.62             	   # fraction of proteins in the organic substrate S1
f_xc    = 0.06                 # fraction of composite in the organic substrate S1
f_xc_pr = 0.2                  # fraction of proteins in the composite Xc
N_S1    = (f_aa+f_pr+f_xc*f_xc_pr)*N_aa

# Influent

pH_in = 6.8

S1_in = T2["S1in"]                             # [g/L]    - COD in solution (soluble COD): S1
S2_in = T2["S2in"]                            # [mmol/L] - Volatile Fatty Acids: S2
C_in  = T2["Cin"]                              # [mmol/L] - Inorganic Carbon
XT_in = T2["XTin"]                              # [g/L]    - Particulate COD in solution : XT
N_in  = T2["Nin"]                              # [mmol/L] - Inorganic Nitrogen

B_in  = Kb*C_in/(Kb+10**(-pH_in))
Z_in  = B_in + S2_in*Ka/(Ka+10**(-pH_in))+N_in

y_in  = [S1_in, S2_in, C_in, Z_in, XT_in]
