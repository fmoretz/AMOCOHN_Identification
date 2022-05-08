# Model Equations
import numpy as np

# Inorganic Carbon evaluation for data comparison
def Ceval(D, X1, X2, qC, Cin, alpha, k4, k5) -> float:
    return  k4*(alpha*D)*X1/D + k5*(alpha*D)*X2/D - qC/D + Cin

# Carbon Dioxide Flow evaluation for data comparison
def qCeval(kLa, C, pH, kH, Pc, Kb) -> float:
    pKb = - np.log10(Kb)
    if type(pH) is not list:
        exp =  float(pH) - float(pKb)
        return kLa*(C/(1+10**exp) - kH*Pc)
    
    else:
        exp = np.empty(shape=(len(pH),))
        ret = np.empty(shape=(len(pH),))
        for i in range(0, len(pH)):
            exp[i] = float(pH[i]) - float(pKb)
            ret[i] = kLa*(C[i]/(1+10**exp[i]) - kH*Pc[i])
        return ret
    
# Particulate evaluation for data comparison
def XTeval(khyd, D, XTin) -> float:
    return D/(D+khyd)*XTin

# Methane Flow evaluation for data comparison
def qMeval(k6, alpha, D, X2) -> float:
    return k6*alpha*D*X2

# Substrate 1 evaluation for data comparison
def S1eval(k1, alpha, D, X1, khyd, Sin, XT) -> float:
    return Sin - k1*alpha*D*X1/D + khyd*XT/D

# Substrate 2 evaluation for data comparison
def S2eval(k2, k3, alpha, D, X2, X1, Sin) -> float:
    return Sin + k2*alpha*D*X1/D - k3*alpha*D*X2/D