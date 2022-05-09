# define the function for RSE
import numpy as np
import solver.RegModel as rm
from typing import List
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error

def RMSE(Data: List, Model: List) -> float:
    MSE = mean_squared_error(Data, Model)
    RMSE = (MSE)**0.5
    return RMSE

def Methane_RMSE(k6, alpha, D, X2, qM):
    """
    Equation for methane flow
    qM / X2 = k6 * alpha * D
    y = m * x
    y = qM / X2 | x = alpha * D | m = k6
    """
    
    y  = np.empty(len(X2))
    m  = np.empty(len(X2))
    x  = np.empty(len(X2))
    mx = np.empty(len(X2))
    
    for i in range(0, len(X2)):
        y[i] = qM[i] / X2[i]
        m[i] = k6 * alpha
        x[i] = D[i]
        mx[i] = m[i] * x[i]
    
    return RMSE(y, mx)

def Methane_LIN(alpha, D, X2, qM):
    """
    Equation for methane flow
    qM / X2 = k6 * alpha * D
    y = m * x
    y = qM / X2 | x = alpha * D | m = k6
    """
    
    Y = np.empty(len(X2))
    X = np.empty(len(X2))
    
    for i in range(0, len(X2)):
        Y[i] = qM[i] / X2[i]
        X[i] = alpha * D[i]
    
    [slope, intercept, r, p, se] = linregress(X, Y)
    
    return [slope, intercept, r**0.5, p, se]

def Carbon_Dioxide_RMSE(kLa, C, CO2, KH, PC, qC):
    """ 
    Equation for carbon dioxide flow
    qC = kLa * (CO2 - KH * PC)
    y = m * x
    y = qC | x = CO2 - KH * Pc | m = kLa
    """
    
    y  = np.empty(len(C))
    m  = np.empty(len(C))
    x  = np.empty(len(C))
    mx = np.empty(len(C))
    
    for i in range(0, len(C)):
        y[i]  = qC[i]
        m[i]  = kLa
        x[i]  = CO2[i]  - KH * PC[i]
        mx[i] = m[i] * x[i]
        
    
    return RMSE(y, mx)

def Carbon_Dioxide_LIN(CO2, PC, KH, qC):
    """ 
    Equation for carbon dioxide flow
    qC = kLa * (CO2 - KH * PC)
    y = m * x
    y = qC | x = CO2 - KH * Pc | m = kLa
    """
    
    X = np.empty(len(PC))
    Y = np.empty(len(PC))

    for i in range(0, len(PC)):
        X[i] = CO2[i] - KH*PC[i]
        Y[i] = qC[i]
             
    [slope, intercept, r, p, se] = linregress(X, Y)
    return [slope, intercept, r**0.5, p, se] 


def Hydrolysis_RMSE(khyd, D, XTin, XT):
    """
    Equation for hydrolysis
    D*(XTin - XT) = khyd * XT
    y = m * x
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

def Hydrolysis_LIN(D, XTin, XT):
    """
    Equation for hydrolysis
    D*(XTin - XT) = khyd * XT
    y = m * x
    y = D*(XTin - XT) | x = XT | m = khyd
    """
    
    X = np.empty(len(XT))
    Y = np.empty(len(XT))

    for i in range(0, len(XT)):
        Y[i] = D[i] * (XTin - XT[i])
        X[i] = XT[i]
             
    [slope, intercept, r, p, se] = linregress(X, Y)
    return [slope, intercept, r**0.5, p, se] 

def Acidogenesis_RMSE(k1, alpha, D, S1, X1, khyd, Sin, XT):
    """
    Equation for acidogenesis
    D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    y = m * x
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
    y = m * x
    y = D*(Sin - S1) + khyd * XT | x = alpha * D * X1 | m = k1
    """
    
    Y = np.empty(len(S1))
    X = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        Y[i] = D[i]*(Sin - S1[i]) + khyd * XT[i]
        X[i] = alpha * D[i] * X1[i]
    
    [slope, intercept, r, p, se] = linregress(X, Y)
    
    return [slope, intercept, r**0.5, p, se]