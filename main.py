# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
from sklearn.metrics import r2_score

# import data and define the variables
from solver.utility.dataimport import *
from solver.utility.PhysConstants import *
from solver.utility.ReactorConf import *
from solver.utility.PostProcess import *
import solver.eqRMSE as solver
import solver.RegModel as rm

# create main function
def main():
    
    # define the bounds
    bd1 = Bounds([0], [5000])
    bd2 = Bounds([0, 0], [5000, 5000])
    bd3 = Bounds([0, 0, 0], [5000, 5000, 5000]) 
    # define the initial values
    v1 = (1,)
    v2 = (1, 2,)
    v3 = (1, 2, 3,)
    
    # define the contraints 
    r2_CH4 = {'type': 'ineq','fun': lambda y: r2_score(qM, rm.qMeval(y, alpha, D, X2))}
    r2_XT  = {'type': 'ineq','fun': lambda y: r2_score(XT, rm.XTeval(y, D, XTin))}
    r2_CO2 = {'type': 'ineq','fun': lambda y: r2_score(qC, rm.qCeval(y, C, pH, KH, PC, Kb))}
    r2_AC  = {'type': 'ineq','fun': lambda y: r2_score(S1, rm.ACeval(y, alpha, D, X1, XT_min.x, S1in, XT))}
    
    # Single Parameter Regression+
    # 'RMSE' Regression
    CH4_min = minimize(solver.Methane_RMSE,    v1, args = (alpha, D, X2, qM), bounds = bd1, tol = 1e-10, constraints = r2_CH4)
    XT_min  = minimize(solver.Hydrolysis_RMSE, v1, args = (D, XTin, XT),      bounds = bd1, tol = 1e-10, constraints = r2_XT)
    
    k6 = CH4_min.x
    khyd = XT_min.x

    # 'LIN' Regression
    [k1, q_AC, score_AC]    = solver.Acidogenesis_LIN(alpha, D, S1, X1, khyd, S1in, XT)
    [kLa, q_CO2, score_CO2] = solver.Carbon_Dioxide_LIN(CO2, PC, KH, qC)
    [[a, b], q_1, score_1]  = solver.Substrate1_LIN(S1, alpha, D)
    [[c, d, e, f], q_2, score_2]  = solver.Substrate2_LIN(S2, alpha, D)
    
    # Parameter retrieving
    mu1m = 1/a
    Ks1  = q_1/0.11
    
    
    mu2m = 1/c
    Ks2  = q_2/0.11
    KI   = 1/(e*mu2m)
    
    # check_11 = mu1m - Ks1/b
    # check_12 = Ks1 - mu1m*b
    # check_21 = KI - 0.11/f
    # check_22 = Ks2 - mu2m*d
    
    # Model evaluation
    qM_Model = np.empty(len(D))
    qC_Model = np.empty(len(D))
    XT_Model = np.empty(len(D))
    AC_Model = np.empty(len(D))
    
    for i in range(0, len(D)):
        qM_Model[i] = rm.qMeval(k6, alpha, D[i], X2[i])
        qC_Model[i] = rm.qCeval(kLa, C[i], pH[i], KH, PC[i], Kb)
        XT_Model[i] = rm.XTeval(khyd, D[i], XTin)
        AC_Model[i] = rm.S1eval(k1, alpha, D[i], X1[i], khyd, S1in, XT[i])

    # print out results
    fprint('RMSE', 'k6', k6, flag = CH4_min.success, r2 = r2_CH4['fun'](k6), itr=CH4_min.nit, funcall=solver.RMSE(qM, rm.qMeval(k6, alpha, D, X2)))
    fprint('RMSE', 'khyd', khyd, flag = XT_min.success, r2 = r2_XT['fun'](khyd), itr=XT_min.nit, funcall=solver.RMSE(XT, rm.XTeval(khyd, D, XTin)))
    fprint('LIN', 'k1', k1, score_AC)
    fprint('LIN', 'kLa', kLa, score_CO2)
    fprint('LIN', 'mu1m + Ks1', [mu1m, Ks1], score_1)
    fprint('LIN', 'mu2m + Ks2 + KI', [mu2m, Ks2, KI], score_2)
    
    # Visualization
    Y_data_CH4  = np.empty(len(D))
    Y_data_CO2  = np.empty(len(D))
    Y_data_AC   = np.empty(len(D))
    Y_data_XT   = np.empty(len(D))
    Y_data_S1   = np.empty(len(D))
    Y_data_S2   = np.empty(len(D))
    
    Y_model_CH4 = np.empty(len(D))
    Y_model_CO2 = np.empty(len(D))
    Y_model_AC  = np.empty(len(D))
    Y_model_XT  = np.empty(len(D))
    Y_model_S1  = np.empty(len(D))
    Y_model_S2  = np.empty(len(D))
    
    X_CH4  = np.empty(len(D))
    X_CO2  = np.empty(len(D))
    X_AC   = np.empty(len(D))
    X_XT   = np.empty(len(D))
    
    X_S1_1   = np.empty(len(D))
    X_S1_2   = np.empty(len(D))
    X_S2_1   = np.empty(len(D))
    X_S2_2   = np.empty(len(D))
    X_S2_3   = np.empty(len(D))
    X_S2_4   = np.empty(len(D))
        
    for i in range(0, len(D)):
        # Methane Flowrate Model - qM / X2 = k6 * alpha * D
        X_CH4[i]        = alpha * D[i]
        Y_model_CH4[i]  = k6 * X_CH4[i]
        Y_data_CH4[i]   = qM[i] / X2[i]
        
        # Carbon Dioxide Flowrate - qC = kLa * (CO2 - KH * PC)
        X_CO2[i]        = CO2[i] - KH*PC[i] 
        Y_model_CO2[i]  = kLa * X_CO2[i]
        Y_data_CO2[i]   = qC[i]
        
        # Hydrolysis - D*(XTin - XT) = khyd * XT
        X_XT[i]         = XT[i]
        Y_model_XT[i]   = khyd * X_XT[i]
        Y_data_XT[i]    = D[i] * (XTin - XT[i])
        
        # Acidogenesis - D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
        X_AC[i]         = alpha * D[i] * X1[i]
        Y_data_AC[i]    = D[i]*(S1in - S1[i]) + khyd * XT[i] 
        Y_model_AC[i]   = k1 * X_AC[i]

        # Susbstrate 1 - S1 = 1/mu1m * alpha * D * S1 / 0.9 + Ks1/mu1m * alpha * D / 0.9 + 0.11 * Ks1
        X_S1_1[i]       = alpha * D[i] * S1[i] / 0.9
        X_S1_2[i]       = alpha * D[i] / 0.9
        Y_data_S1[i]    = S1[i]
        Y_model_S1[i]   = 1/mu1m * X_S1_1[i] + Ks1/mu1m * X_S1_2[i] + 0.11 * Ks1
        
        # Substrate 2 - S2 = alpha / (0.9 * mu2m) * D * (S2 + Ks2 + S2^2/KI) + 0.11 * (Ks2 + S2^2/KI)
        X_S2_1[i]       = alpha * D[i] * S2[i] / 0.9
        X_S2_2[i]       = alpha * D[i] / 0.9
        X_S2_3[i]       = alpha * D[i] * S2[i]**2 / 0.9
        X_S2_4[i]       = S2[i]**2
        Y_data_S2[i]    = S2[i]
        Y_model_S2[i]   = 1/mu2m * X_S2_1[i] + Ks2/mu2m * X_S2_2[i] + 1/( mu2m * KI ) * X_S2_3[i] + 0.11 / KI * X_S2_4[i] + 0.11 * Ks2

    regplot(X_CH4,  [Y_data_CH4, Y_model_CH4], 'D [1/d]', 'qM/X2 [1/d]', 'k6 regression - Methane Flowrate')
    regplot(X_CO2,  [Y_data_CO2, Y_model_CO2], 'CO2-KH*PC [mol/m3]', 'qC [mol/m3/d]', 'kLa regression - Carbon Dioxide Flowrate')
    regplot(X_XT,   [Y_data_XT,  Y_model_XT],  'XT [kg/m3]', 'D*(XTin-XT) [kg/m3/d]', 'khyd regression - Hydrolysis')
    regplot(X_AC,   [Y_data_AC,  Y_model_AC],  'alpha*D*X1 [kg/m3/d]', 'D*(S1in-S1)+khyd*XT [kg/m3/d]', 'k1 regression - Acidogenesis')
    regplot(X_S1_1, [Y_data_S1, Y_model_S1],   'alpha*D*S1/0.9 [kg/m3/d]', 'S1 [kg/m3]', 'mu1m regression - Substrate 1 - X: 1' )
    regplot(X_S1_2, [Y_data_S1, Y_model_S1],   'alpha*D/0.9 [d-1]', 'S1 [kg/m3]', 'mu1m regression - Substrate 1 - X: 2' )
    
    multiplot([Y_data_S1, Y_model_S1], [X_S1_1, X_S1_2], ['alpha*D*S1/0.9 [kg/m3/d]','alpha*D/0.9 [d-1]'], 'S1 [kg/m3]', 'mu1m and Ks1 regression - Substrate 1', showbool=True)

    return 0;

if __name__ == "__main__":
    main()
    
