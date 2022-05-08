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
    [kLa, q, r2, p_val, err] = solver.Carbon_Dioxide_LIN(C, pH, KH, PC, Kb, qC)
    [k1, q1, r21, p_val1, err1] = solver.Acidogenesis_LIN(alpha, D, S1, X1, XT_min.x, S1in, XT)
    
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
    fprint('k6-RMSE', CH4_min.x, CH4_min.success, CH4_min.nit, CH4_min.fun, [qM, qM_Model])
    fprint('kLa-RMSE', CO2_min.x, CO2_min.success, CO2_min.nit, CO2_min.fun, [qC, qC_Model])
    fprint('kLa-LIN', kLa, 'NO-INFO', err, r2, [qC, rm.qCeval(kLa, C, pH, KH, PC, Kb)])
    fprint('khyd-RMSE', XT_min.x, XT_min.success, XT_min.nit, XT_min.fun, [XT, XT_Model])
    fprint('k1-RMSE', S1_min.x, S1_min.success, S1_min.nit, S1_min.fun, [S1, S1_Model])
    fprint('k1-LIN', k1, 'NO-INFO', err, r2, [S1, S1_Model])

    # plot results
    legend = ['Experimental data', 'Model']
    regplot(1/D, [qM, qM_Model], 'time [days]', 'qM [g/s]', 'Methane production regression', legend, False)
    regplot(1/D, [qC, qC_Model], 'time [days]', 'qC [g/s]', 'Carbon Dioxide production regression', legend, False)
    regplot(1/D, [XT, XT_Model], 'time [days]', 'XT [kg/m3]', 'Particulate Hydrolysis regression', legend, False)
    regplot(1/D, [S1, S1_Model], 'time [days]', 'S1 [kg/m3]', 'Acidogenesis regression', legend, True)

    return 0;

if __name__ == "__main__":
    main()
    
