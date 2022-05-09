# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
from sklearn.metrics import r2_score
plt.close('all')

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
    v1 = 20
    v2 = [20, 20]
    v3 = [20, 20, 20]
    
    # define the contraints 
    r2_CH4 = {'type': 'ineq','fun': lambda y: r2_score(qM, rm.qMeval(y, alpha, D, X2))}
    r2_XT  = {'type': 'ineq','fun': lambda y: r2_score(XT, rm.XTeval(y, D, XTin))}
    r2_CO2 = {'type': 'ineq','fun': lambda y: r2_score(qC, rm.qCeval(y, C, pH, KH, PC, Kb))}
    r2_S1  = {'type': 'ineq','fun': lambda y: r2_score(S1, rm.S1eval(y, alpha, D, X1, XT_min.x, S1in, XT))}
    
    # Single Parameter Regression+
    # 'RMSE' Regression
    CH4_min = minimize(solver.Methane_RMSE,    x0 = (v1,), args = (alpha, D, X2, qM), bounds = bd1, tol = 1e-10, constraints = r2_CH4)
    XT_min  = minimize(solver.Hydrolysis_RMSE, x0 = (v1,), args = (D, XTin, XT),      bounds = bd1, tol = 1e-10, constraints = r2_XT)
    S1_min  = minimize(solver.Substrate1_RMSE, x0 = (v2,), args = (S1, alpha, D),     bounds = bd2, tol = 1e-10)
    S2_min  = minimize(solver.Substrate2_RMSE, x0 = (v3,), args = (S2, alpha, D),     bounds = bd3, tol = 1e-10)
    
    k6 = CH4_min.x
    khyd = XT_min.x
    [mu1m, Ks1] = S1_min.x
    [mu2m, Ks2, KI] = S2_min.x

    # 'LIN' Regression
    [k1, q_S1, score_S1]      = solver.Acidogenesis_LIN(alpha, D, S1, X1, khyd, S1in, XT)
    [kLa, q_CO2, score_CO2]   = solver.Carbon_Dioxide_LIN(CO2, PC, KH, qC)
    
    
    # Model evaluation
    qM_Model = np.empty(len(D))
    qC_Model = np.empty(len(D))
    XT_Model = np.empty(len(D))
    S1_Model = np.empty(len(D))
    
    for i in range(0, len(D)):
        qM_Model[i] = rm.qMeval(k6, alpha, D[i], X2[i])
        qC_Model[i] = rm.qCeval(kLa, C[i], pH[i], KH, PC[i], Kb)
        XT_Model[i] = rm.XTeval(khyd, D[i], XTin)
        S1_Model[i] = rm.S1eval(k1, alpha, D[i], X1[i], khyd, S1in, XT[i])

    # print out results
    fprint('RMSE', 'k6', k6, flag = CH4_min.success, r2 = r2_CH4['fun'](k6), itr=CH4_min.nit, funcall=solver.RMSE(qM, rm.qMeval(k6, alpha, D, X2)))
    fprint('RMSE', 'khyd', khyd, flag = XT_min.success, r2 = r2_XT['fun'](khyd), itr=XT_min.nit, funcall=solver.RMSE(XT, rm.XTeval(khyd, D, XTin)))
    fprintf('RMSE', 'mu1_max + Ks1', [mu1m, Ks1], flag = S1_min.success, itr=S1_min.nit)
    fprintf('RMSE', 'mu2_max + Ks2 + KI', [mu2m, Ks2, KI], flag = S2_min.success, itr=S2_min.nit)
    fprint('LIN', 'k1', k1, score_S1)
    fprint('LIN', 'kLa', kLa, score_CO2)
    
    # Visualization
    Y_data_CH4  = np.empty(len(D))
    Y_data_CO2  = np.empty(len(D))
    Y_data_S1   = np.empty(len(D))
    Y_data_XT   = np.empty(len(D))
    
    Y_model_CH4 = np.empty(len(D))
    Y_model_CO2 = np.empty(len(D))
    Y_model_S1  = np.empty(len(D))
    Y_model_XT  = np.empty(len(D))
    
    X_CH4  = np.empty(len(D))
    X_CO2  = np.empty(len(D))
    X_S1   = np.empty(len(D))
    X_XT   = np.empty(len(D))
    
    for i in range(0, len(D)):
        # Methane Flowrate Model - qM / X2 = k6 * alpha * D
        X_CH4[i]       = alpha * D[i]
        Y_model_CH4[i] = k6 * X_CH4[i]
        Y_data_CH4[i]  = qM[i] / X2[i]
        
        # Carbon Dioxide Flowrate - qC = kLa * (CO2 - KH * PC)
        X_CO2[i]       = CO2[i] - KH*PC[i] 
        Y_model_CO2[i] = kLa * X_CO2[i]
        Y_data_CO2[i]  = qC[i]
        
        # Hydrolysis - D*(XTin - XT) = khyd * XT
        X_XT[i]        = XT[i]
        Y_model_XT[i]  = khyd * X_XT[i]
        Y_data_XT[i]   = D[i] * (XTin - XT[i])
        
        # Acidogenesis - D*(Sin - S1) + khyd * XT = k1 * alpha * D * X1
        X_S1[i]        = alpha * D[i] * X1[i]
        Y_data_S1[i]   = D[i]*(S1in - S1[i]) + khyd * XT[i] 
        Y_model_S1[i]  = k1 * X_S1[i]

    regplot(X_CH4, [Y_data_CH4, Y_model_CH4], 'D [1/d]', 'qM/X2 [1/d]', 'k6 regression - Methane Flowrate')
    regplot(X_CO2, [Y_data_CO2, Y_model_CO2], 'CO2-KH*PC [mol/m3]', 'qC [mol/m3/d]', 'kLa regression - Carbon Dioxide Flowrate')
    regplot(X_XT,  [Y_data_XT,  Y_model_XT],  'XT [kg/m3]', 'D*(XTin-XT) [kg/m3/d]', 'khyd regression - Hydrolysis')
    regplot(X_S1,  [Y_data_S1,  Y_model_S1],  'alpha*D*X1 [kg/m3/d]', 'D*(S1in-S1)+khyd*XT [kg/m3/d]', 'k1 regression - Acidogenesis', showbool=True)

    return 0;

if __name__ == "__main__":
    main()
    
