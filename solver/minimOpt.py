# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
from sklearn.metrics import r2_score
plt.close('all')

# import data and define the variables
from utility.dataimport import *
from utility.PhysConstants import *
from utility.ReactorConf import *
from utility.PostProcess import *
import eqRMSE as solver
import RegModel as rm

# create main function
def main():
    
    # define the bounds
    bd = Bounds([0], [5000])
    
    # define the initial values
    gs = 20
    
    # define the contraints 
    r2_CH4 = {'type': 'ineq','fun': lambda y: r2_score(qM, rm.qMeval(y, alpha, D, X2))}
    r2_CO2 = {'type': 'ineq','fun': lambda y: r2_score(qC, rm.qCeval(y,  C, pH, KH, PC, Kb))}
    r2_XT  = {'type': 'ineq','fun': lambda y: r2_score(XT, rm.XTeval(y, D, XTin))}
    r2_S1  = {'type': 'ineq','fun': lambda y: r2_score(S1, rm.S1eval(y, alpha, D, X1, XT_min.x, S1in, XT))}
    
    # Single Parameter Regression
    CH4_min = minimize(solver.Methane_RMSE, x0 = (gs,), args = (alpha, D, X2, qM), bounds = bd, tol = 1e-10, constraints=r2_CH4)
    CO2_min = minimize(solver.Carbon_Dioxide_RMSE, x0 = (gs,), args = (C, pH, KH, PC, Kb, qC), bounds = bd, tol = 1e-10, constraints=r2_CO2)
    XT_min  = minimize(solver.Hydrolysis_RMSE, x0 = (gs,), args = (D, XTin, XT), bounds = bd, tol = 1e-10, constraints=r2_XT)
    S1_min  = minimize(solver.Acidogenesis_RMSE, x0 = (gs,), args = (alpha, D, S1, X1, XT_min.x, S1in, XT), bounds = bd, tol = 1e-10, constraints=r2_S1)
    [kLa, q_CO2, r2_CO2, p_CO2, err_CO2] = solver.Carbon_Dioxide_LIN(C, pH, KH, PC, Kb, qC)
    [k1, q_S1, r2_S1, p_S1, err_S1] = solver.Acidogenesis_LIN(alpha, D, S1, X1, XT_min.x, S1in, XT)
    
    # Model evaluation
    qM_Model = np.empty(len(D))
    qC_Model = np.empty(len(D))
    XT_Model = np.empty(len(D))
    S1_Model = np.empty(len(D))
    
    for i in range(0, len(D)):
        qM_Model[i] = rm.qMeval(CH4_min.x, alpha, D[i], X2[i])
        qC_Model[i] = rm.qCeval(kLa, C[i], pH[i], KH, PC[i], Kb)
        XT_Model[i] = rm.XTeval(XT_min.x, D[i], XTin)
        S1_Model[i] = rm.S1eval(k1, alpha, D[i], X1[i], XT_min.x, S1in, XT[i])

    # print out results
    fprint('RMSE', 'k6', CH4_min.x, flag = CH4_min.success, r2 = r2_CH4['fun'](CH4_min.x), itr=CH4_min.nit, funcall=solver.RMSE(qM, rm.qMeval(CH4_min.x, alpha, D, X2)))
    fprint('RMSE', 'khyd', XT_min.x, flag = XT_min.success, r2 = r2_XT['fun'](XT_min.x), itr=XT_min.nit, funcall=solver.RMSE(XT, rm.XTeval(XT_min.x, D, XTin)))
    fprint('LIN', 'kLa', kLa, r2_CO2, err_CO2)
    fprint('LIN', 'k1', k1, r2_S1, err_S1)
    

    return 0;

if __name__ == "__main__":
    main()
    
