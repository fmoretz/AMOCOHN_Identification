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

# model equations
def CH4model(alpha, D, X2, qM, k6):
    """
     Methane Flowrate Model - qM / X2 = k6 * alpha * D
    """
    Y_data_CH4  = np.empty(len(D))
    Y_model_CH4 = np.empty(len(D))
    X_CH4  = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_CH4[i]        = alpha * D[i]
        Y_model_CH4[i]  = k6 * X_CH4[i]
        Y_data_CH4[i]   = qM[i] / X2[i]
        
    return X_CH4, Y_data_CH4, Y_model_CH4
    
def ACmodel(alpha, D, S1in, S1, XT, X1, k1, khyd):
    """
    Acidogenesis - D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    """
    Y_data_AC   = np.empty(len(D))
    Y_model_AC  = np.empty(len(D))
    X_AC   = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_AC[i]         = alpha * D[i] * X1[i]
        Y_data_AC[i]    = D[i]*(S1in - S1[i]) + khyd * XT[i] 
        Y_model_AC[i]   = k1 * X_AC[i]

    return X_AC, Y_data_AC, Y_model_AC

def CO2model(qC, CO2, KH, PC, kLa):
    """
    Carbon Dioxide Flowrate - qC = kLa * (CO2 - KH * PC)
    """
    Y_data_CO2  = np.empty(len(qC))
    Y_model_CO2 = np.empty(len(qC))
    X_CO2  = np.empty(len(qC))
    
    for i in range(0, len(qC)):
        X_CO2[i]        = CO2[i] - KH*PC[i] 
        Y_model_CO2[i]  = kLa * X_CO2[i]
        Y_data_CO2[i]   = qC[i]
    
    return X_CO2, Y_data_CO2, Y_model_CO2

def XTmodel(D, XTin, XT, khyd):
    """
    Hydrolysis - D*(XTin - XT) = khyd * XT
    """
    Y_data_XT   = np.empty(len(D))
    Y_model_XT  = np.empty(len(D))
    X_XT   = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_XT[i]         = XT[i]
        Y_model_XT[i]   = khyd * X_XT[i]
        Y_data_XT[i]    = D[i] * (XTin - XT[i])
        
    return X_XT, Y_data_XT, Y_model_XT
    
    
def S1model(alpha, D, mu1m, Ks1, S1):
    """
    Susbstrate 1 - S1 = 1/mu1m * alpha * D * S1 / 0.9 + Ks1/mu1m * alpha * D / 0.9 + 0.11 * Ks1
    """
    Y_data_S1   = np.empty(len(D))
    Y_model_S1  = np.empty(len(D))
    X_S1_1   = np.empty(len(D))
    X_S1_2   = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_S1_1[i]       = alpha * D[i] * S1[i] / 0.9
        X_S1_2[i]       = alpha * D[i] / 0.9
        Y_data_S1[i]    = S1[i]
        Y_model_S1[i]   = 1/mu1m * X_S1_1[i] + Ks1/mu1m * X_S1_2[i] + 0.11 * Ks1
        
    return X_S1_1, X_S1_2, Y_data_S1, Y_model_S1
    
    
def S2model(alpha, D, mu2m, Ks2, KI, S2):
    """
    Substrate 2 - S2 = alpha / (0.9 * mu2m) * D * (S2 + Ks2 + S2^2/KI) + 0.11 * (Ks2 + S2^2/KI)
    """
    Y_model_S2  = np.empty(len(D))    
    Y_data_S2   = np.empty(len(D))
    X_S2_1   = np.empty(len(D))
    X_S2_2   = np.empty(len(D))
    X_S2_3   = np.empty(len(D))
    X_S2_4   = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_S2_1[i]       = alpha * D[i] * S2[i] / 0.9
        X_S2_2[i]       = alpha * D[i] / 0.9
        X_S2_3[i]       = alpha * D[i] * S2[i]**2 / 0.9
        X_S2_4[i]       = S2[i]**2
        Y_data_S2[i]    = S2[i]
        Y_model_S2[i]   = 1/mu2m * X_S2_1[i] + Ks2/mu2m * X_S2_2[i] + 1/( mu2m * KI ) * X_S2_3[i] + 0.11 / KI * X_S2_4[i] + 0.11 * Ks2

    return X_S2_1, X_S2_2, X_S2_3, X_S2_4, Y_data_S2, Y_model_S2
    
def MEmodel(alpha, D, X1, X2, S2, S2in, k2, k3):
    """
    Methanogenesis - D*(S2in - S2) = k3 * alpha * D * X2 - k2 * alpha * D * X1
    """
    Y_data_ME   = np.empty(len(D))
    Y_model_ME  = np.empty(len(D))
    X_ME_1      = np.empty(len(D))
    X_ME_2      = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_ME_1[i]       = alpha * D[i] * X1[i]
        X_ME_2[i]       = alpha * D[i] * X2[i]
        Y_data_ME[i]    = D[i]*(S2in - S2[i]) 
        Y_model_ME[i]   = k3 * X_ME_2[i] - k2 * X_ME_1[i]

    return X_ME_1, X_ME_2, Y_data_ME, Y_model_ME

def ICmodel(alpha, D, qC, X1, X2, C, Cin, k5, k4):
    """
    Carbogenesis - qC - D*(Cin - C) = k5 * alpha * D * X2 - k4 * alpha * D * X1
    """
    Y_data_IC   = np.empty(len(D))
    Y_model_IC  = np.empty(len(D))
    X_IC_1      = np.empty(len(D))
    X_IC_2      = np.empty(len(D))
    
    for i in range(0, len(D)):
        X_IC_1[i]       = alpha * D[i] * X1[i]
        X_IC_2[i]       = alpha * D[i] * X2[i]
        Y_data_IC[i]    = qC[i] - D[i]*(Cin - C[i]) 
        Y_model_IC[i]   = k5 * X_IC_2[i] - k4 * X_IC_1[i]

    return X_IC_1, X_IC_2, Y_data_IC, Y_model_IC