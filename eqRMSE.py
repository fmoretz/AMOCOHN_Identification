# define the function for RSE
import numpy as np
from typing import List
import RegModel as rm
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error

def RMSE(Data: List, Model: List) -> float:
    MSE = mean_squared_error(Data, Model)
    RMSE = (MSE)**0.5
    return RMSE

def Methane_RMSE(k6, alpha, D, X2, qM):
    """
    Equation for methane flow
    qM = k6 * alpha * D * X2
    y = m * x + b
    y = qM | x = X2 | m = k6 * alpha * D
    k6 = m / alpha * D
    """
    
    y  = np.empty(len(X2))
    m  = np.empty(len(X2))
    x  = np.empty(len(X2))
    mx = np.empty(len(X2))
    
    for i in range(0, len(X2)):
        y[i] = qM[i]
        m[i] = k6 * alpha * D[i]
        x[i] = X2[i]
        mx[i] = m[i] * x[i]
    
    return RMSE(y, mx)

def Carbon_Dioxide_RMSE(kLa, C, pH, KH, PC, Kb, qC):
    """ 
    Equation for carbon dioxide flow
    qC = kLa * (C / (1 + 10**(pH - pKb)) - kH * Pc)
    y = m * x + b
    y = qC | x = C / (1 + 10**(pH - pKb)) - kH * Pc | m = kLa
    """
    
    y  = np.empty(len(C))
    m  = np.empty(len(C))
    x  = np.empty(len(C))
    mx = np.empty(len(C))
    
    pKb = -np.log10(Kb)
    
    for i in range(0, len(C)):
        y[i]  = qC[i]
        m[i]  = kLa
        x[i]  = C[i] * (1 / (1 + 10**(pH[i] - pKb))) - KH * PC[i]
        mx[i] = m[i] * x[i]
    
    return RMSE(y, mx)

def Carbon_Dioxide_LIN(C, pH, KH, PC, Kb, qC):
    """ 
    Equation for carbon dioxide flow
    qC = kLa * (C / (1 + 10**(pH - pKb)) - kH * Pc)
    y = m * x + b
    y = qC | x = C / (1 + 10**(pH - pKb)) - kH * Pc | m = kLa
    """
    pKb = -np.log10(Kb)
    Y = qC
    
    if type(pH) is not list:
        exp =  float(pH) - float(pKb)
        X = C*1/(1+10**(pH-pKb)) - KH*PC
    
    else:
        exp = np.empty(shape=(len(pH),))
        X = np.empty(shape=(len(pH),))
        for i in range(0, len(pH)):
            exp[i] = float(pH[i]) - float(pKb)
            X[i] = (C[i]/(1+10**exp[i]) - KH*PC[i])

    [slope, intercept, r, p, se] = linregress(X, Y)
    return [slope, intercept, r**0.5, p, se] 


def Hydrolysis_RMSE(khyd, D, XTin, XT):
    """
    Equation for hydrolysis
    D*(XTin - XT) = khyd * XT
    y = m * x + b
    y = D*(XTin - XT) | x = XT | m = khyd
    """
    
    y = np.empty(len(XT))
    m = np.empty(len(XT))
    x = np.empty(len(XT))
    mx = np.empty(len(XT))
    
    for i in range(0, len(XT)):
        y[i]  = D[i] * (XTin - XT[i])
        x[i]  = XT[i]
        m[i]  = khyd
        mx[i] = m[i] * x[i]

    return RMSE(y, mx)

def Acidogenesis_RMSE(k1, alpha, D, S1, X1, khyd, Sin, XT):
    """
    Equation for acidogenesis
    D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    y = m * x + b
    y = D*(Sin - S1) + khyd * XT | x = alpha * D * X1 | m = k1
    """
    
    y = np.empty(len(S1))
    x = np.empty(len(S1))
    m = np.empty(len(S1))
    mx = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        y[i] = D[i] * (Sin - S1[i]) + khyd * XT[i]
        x[i] = alpha * D[i] * X1[i]
        m[i] = k1
        mx[i] = m[i] * x[i]
    
    return RMSE(y, mx)

def Acidogenesis_LIN(alpha, D, S1, X1, khyd, Sin, XT):
    """
    Equation for acidogenesis
    D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    y = m * x + b
    y = D*(Sin - S1) + khyd * XT | x = alpha * D * X1 | m = k1
    """
    
    Y = np.empty(len(S1))
    X = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        Y[i] = D[i]*(Sin - S1[i]) + khyd * XT[i]
        X[i] = alpha * D[i] * X1[i]
    
    [slope, intercept, r, p, se] = linregress(X, Y)
    
    return [slope, intercept, r**0.5, p, se]