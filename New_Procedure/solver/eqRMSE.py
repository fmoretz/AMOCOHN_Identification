# define the function for RSE
import numpy as np
from typing import List
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
import pandas as pd

def RMSE(Data: List, Model: List) -> float:
    """
    RMSE tool for regression model
    RMSE = sqrt(sum((Data - Model)^2) / len(Data))
    """
    
    MSE = mean_squared_error(Data, Model)
    RMSE = (MSE)**0.5
    return RMSE

# --------------------------------------------------------------------------------------------------------------------- RMSE Algorithm
def Substrate1_RMSE(x, S1, alpha, D):
    """
    Equation for substrate 1
    S1 = alpha / (0.9 * mu1m) * D * (S1 + Ks1) + 0.11 * Ks1
    y = a * x1 + b * x2 + q
    y = S1 | a = 1 / mu1m | x1 = alpha * D * S1/ 0.9 | 
    b = Ks1 / mu1m | x2 = alpha * D / 0.9 | 
    q = 0.11 * Ks1
    """
    
    mu1m, Ks1 = x
    
    y  = np.empty(len(S1))
    a  = np.empty(len(S1))
    x1 = np.empty(len(S1))
    b  = np.empty(len(S1))
    x2 = np.empty(len(S1))
    q  = np.empty(len(S1))
    ax1bx2q = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        y[i]  = S1[i]
        a[i]  = 1 / mu1m
        x1[i] = alpha * D[i] * S1[i] / 0.9
        b[i]  = Ks1 / mu1m
        x2[i] = alpha * D[i] / 0.9
        q[i]  = 0.11 * Ks1
        ax1bx2q[i] = a[i] * x1[i] + b[i] * x2[i] + q[i]
        
    return RMSE(y, ax1bx2q)

def Substrate2_RMSE(mu2m, Ks2, KI, S2, alpha, D):
    """
    Equation for substrate 2
    S2 = alpha / (0.9 * mu2m) * D * (S2 + Ks2 + S2^2/KI) + 0.11 * (Ks2 + S2^2/KI)
    y = a * x1 + b * x2 + c * x3 + d * x4 + q
    y = S1 | a = 1 / mu2m | x1 = alpha * D * S2/ 0.9 |
    b = Ks2 / mu2m | x2 = alpha * D / 0.9 | c = 1 / (mu2m * KI) |
    x3 = alpha * D * S2^2 / 0.9 | d = 0.11 / KI | x4 = S2^2 | 
    q = 0.11 * Ks2
    """
    y  = np.empty(len(S2))
    a  = np.empty(len(S2))
    x1 = np.empty(len(S2))
    b  = np.empty(len(S2))
    x2 = np.empty(len(S2))
    c  = np.empty(len(S2))
    x3 = np.empty(len(S2))
    d  = np.empty(len(S2))
    x4 = np.empty(len(S2))
    q  = np.empty(len(S2))
    ax1bx2cx3dx4q = np.empty(len(S2))
    
    for i in range(0, len(S2)):
        y[i]  = S2[i]
        a[i]  = 1 / mu2m[i]
        x1[i] = alpha * D[i] * S2[i] / 0.9
        b[i]  = Ks2[i] / mu2m[i]
        x2[i] = alpha * D[i] / 0.9
        c[i]  = 1 / (mu2m[i] * KI[i])
        x3[i] = alpha * D[i] * S2[i] ** 2 / 0.9
        d[i]  = 0.11 / KI[i]
        x4[i] = S2[i] ** 2
        q[i]  = 0.11 * Ks2[i]
        ax1bx2cx3dx4q[i] = a[i] * x1[i] + b[i] * x2[i] + c[i] * x3[i] + d[i] * x4[i] + q[i]

    return RMSE(y, ax1bx2cx3dx4q)

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
        y[i]  = qM[i] / X2[i]
        m[i]  = k6 * alpha
        x[i]  = D[i]
        mx[i] = m[i] * x[i]
    
    return RMSE(y, mx)

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

def Acidogenesis_RMSE(k1, alpha, D, S1, X1, khyd, Sin, XT):
    """
    Equation for acidogenesis
    D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    y = m * x + q
    y = D*(Sin - S1) + khyd * XT | x = alpha * D * X1 | m = k1 | q = khyd
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

# --------------------------------------------------------------------------------------------------------------------- LIN Algorithm

def Methane_LIN(alpha, D, X2, qM):
    """
    Equation for methane flow
    qM / X2 = k6 * alpha * D
    y = m * x
    y = qM / X2 | x = alpha * D | m = k6
    """
    
    Y = np.empty(len(X2)).reshape(-1,1)
    X = np.empty(len(X2)).reshape(-1,1)
    
    for i in range(0, len(X2)):
        Y[i] = qM[i] / X2[i]
        X[i] = alpha * D[i]
        
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Carbon_Dioxide_LIN(CO2, PC, KH, qC):
    """ 
    Equation for carbon dioxide flow
    qC = kLa * (CO2 - KH * PC)
    y = m * x
    y = qC | x = CO2 - KH * Pc | m = kLa
    """
    
    X = np.empty(len(PC)).reshape(-1,1)
    Y = np.empty(len(PC)).reshape(-1,1)

    for i in range(0, len(PC)):
        X[i] = CO2[i] - KH*PC[i]
        Y[i] = qC[i]
             
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Hydrolysis_LIN(D, XTin, XT):
    """
    Equation for hydrolysis
    D*(XTin - XT) = khyd * XT
    y = m * x
    y = D*(XTin - XT) | x = XT | m = khyd
    """
    
    X = np.empty(len(XT)).reshape(-1,1)
    Y = np.empty(len(XT)).reshape(-1,1)

    for i in range(0, len(XT)):
        Y[i] = D[i] * (XTin - XT[i])
        X[i] = XT[i]
             
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Acidogenesis_LIN(alpha, D, S1, X1, khyd, Sin, XT):
    """
    Equation for acidogenesis
    D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
    y = m * x
    y = D*(Sin - S1) + khyd * XT | x = alpha * D * X1 | m = k1
    """
    
    Y = np.empty(len(S1)).reshape(-1,1)
    X = np.empty(len(S1)).reshape(-1,1)
    
    for i in range(0, len(S1)):
        Y[i] = D[i]*(Sin - S1[i]) + khyd * XT[i]
        X[i] = alpha * D[i] * X1[i]
    
    reg_ = LinearRegression(fit_intercept=False).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Substrate1_LIN(S1, alpha, D):
    """
    Equation for substrate 1
    S1 = alpha / (0.9 * mu1m) * D * (S1 + Ks1) + 0.11 * Ks1
    y = a * x1 + b * x2 + q
    y = S1 | a = 1 / mu1m | x1 = alpha * D * S1/ 0.9 | 
    b = Ks1 / mu1m | x2 = alpha * D / 0.9 | 
    q = 0.11 * Ks1
    """
    
    y  = np.empty(len(S1))
    x1 = np.empty(len(S1))
    x2 = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        y[i]  = S1[i]
        x1[i] = alpha * D[i] * S1[i] / 0.9
        x2[i] = alpha * D[i] / 0.9
    
    data = {'x1': x1, 'x2': x2, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'y'])
    Y = df['y']
    X = df[['x1', 'x2']]
        
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Substrate2_LIN(S2, alpha, D):
    """
    Equation for substrate 2
    S2 = alpha / (0.9 * mu2m) * D * (S2 + Ks2 + S2^2/KI) + 0.11 * (Ks2 + S2^2/KI)
    y = a * x1 + b * x2 + c * x3 + d * x4 + q
    y = S1 | a = 1 / mu2m | x1 = alpha * D * S2/ 0.9 |
    b = Ks2 / mu2m | x2 = alpha * D / 0.9 | c = 1 / (mu2m * KI) |
    x3 = alpha * D * S2^2 / 0.9 | d = 0.11 / KI | x4 = S2^2 | 
    q = 0.11 * Ks2
    """
    y  = np.empty(len(S2))
    x1 = np.empty(len(S2))
    x2 = np.empty(len(S2))
    x3 = np.empty(len(S2))
    x4 = np.empty(len(S2))
    
    for i in range(0, len(S2)):
        y[i]  = S2[i]
        x1[i] = alpha * D[i] * S2[i] / 0.9
        x2[i] = alpha * D[i] / 0.9
        x3[i] = alpha * D[i] * S2[i] ** 2 / 0.9
        x4[i] = S2[i] ** 2

    data = {'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'x3', 'x4', 'y'])
    Y = df['y']
    X = df[['x1', 'x2', 'x3', 'x4']]
        
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Methanogenesis_LIN(alpha, D, X1, X2, S2, S2in):
    """
    Equation for methanogenesis
    D*(S2in - S2) = k3 * alpha * D * X2 - k2 * alpha * D * X1
    y = a * x1 + b * x2
    y = D*(S2in - S2) | a = k3 | x1 = alpha * D * X2
    b = k2 | x2 = - alpha * D * X1 
    """
    
    y  = np.empty(len(X1))
    x1 = np.empty(len(X1))
    x2 = np.empty(len(X1))
    
    for i in range(0, len(X1)):
        y[i]  = D[i]*(S2in - S2[i])
        x1[i] = alpha * D[i] * X2[i]
        x2[i] = - alpha * D[i] * X1[i]
    
    data = {'x1': x1, 'x2': x2, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'y'])
    Y = df['y']
    X = df[['x1', 'x2']]
        
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Carbogenesis_LIN(alpha, D, qC, X1, X2, C, Cin):
    """
    Equation for inorganic carbon production (Carbogenesis)
    qC - D*(Cin - C) = k5 * alpha * D * X2 - k4 * alpha * D * X1
    y = a * x1 + b * x2
    y = qC - D*(Cin - C) | a = k5 | x1 = alpha * D * X2
    b = k4 | x2 = - alpha * D * X1 
    """
    
    y   = np.empty(len(X1))
    x1  = np.empty(len(X1))
    x2  = np.empty(len(X1))
    
    for i in range(0, len(X1)):
        y[i]  = qC[i] - D[i]*(Cin - C[i])
        x1[i] = alpha * D[i] * X2[i]
        x2[i] = - alpha * D[i] * X1[i]
        
    data = {'x1': x1, 'x2': x2, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'y'])
    Y = df['y']
    X = df[['x1', 'x2']]
    
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Methanogenesis_RATIO(D, qM, S2in, S2, S1in, S1, khyd, XT):
    """
    New expression for methane flowrate evaluation:
    qM = k6/k3 * D*(S2in - S2) + k6 * k2/(k3 * k1) * (D*(Sin - S1) + khyd * XT)
    y = a * x1 + b * x2
    y = qM | a = k6/k3 | x1 = D*(S2in - S2)
    b = k6 * k2/(k3 * k1) | x2 = (D*(Sin - S1) + khyd * XT)
    """
    
    y   = np.empty(len(S1))
    x1  = np.empty(len(S1))
    x2  = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        y[i]  = qM[i]
        x1[i] = D[i] * (S2in - S2[i])
        x2[i] = D[i] * (S1in - S1[i]) + khyd * XT[i]
        
    data = {'x1': x1, 'x2': x2, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'y'])
    Y = df['y']
    X = df[['x1', 'x2']]
    
    reg_ = LinearRegression(fit_intercept=True).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]

def Carbogenesis_RATIO(D, qC, qM, Cin, C, S1in, S1, khyd, XT):
    """
    New expression for carbon dioxide flowrate evaluation:
    qC - D*(Cin - C) = k4/k1 * (D*(Sin - S1) + khyd * XT) + k5 / k6 * qM
    y = a * x1 + b * x2
    y = qC - D*(Cin - C) | a = k4/k1 | x1 = (D*(Sin - S1) + khyd * XT)
    b = k5 / k6 | x2 = qM
    """
    
    y   = np.empty(len(S1))
    x1  = np.empty(len(S1))
    x2  = np.empty(len(S1))
    
    for i in range(0, len(S1)):
        y[i]  = qC[i] - D[i]*(Cin - C[i])
        x1[i] = D[i] * (S1in - S1[i]) + khyd * XT[i]
        x2[i] = qM[i]
        
    data = {'x1': x1, 'x2': x2, 'y': y}
    df = pd.DataFrame(data, columns=['x1', 'x2', 'y'])
    Y = df['y']
    X = df[['x1', 'x2']]
    
    reg_ = LinearRegression(fit_intercept=False).fit(X, Y)
    return [reg_.coef_, reg_.intercept_, reg_.score(X, Y)]