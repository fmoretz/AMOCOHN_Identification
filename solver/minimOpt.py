# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
from sklearn.metrics import r2_score
from utility.PhysConstants import *
from utility.ReactorConf import *
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
    r2_XT  = {'type': 'ineq','fun': lambda y: r2_score(XT, rm.XTeval(y, D, XTin))}
    
    # Single Parameter Regression
    CH4_min = minimize(solver.Methane_RMSE, x0 = (gs,), args = (alpha, D, X2, qM), bounds = bd, tol = 1e-10, constraints=r2_CH4)
    XT_min  = minimize(solver.Hydrolysis_RMSE, x0 = (gs,), args = (D, XTin, XT), bounds = bd, tol = 1e-10, constraints=r2_XT)
    [kLa, kH, r2_CO2, p_CO2, err_CO2] = solver.Carbon_Dioxide_LIN(C, pH, PC, Kb, qC)
    [k1, q_S1, r2_S1, p_S1, err_S1] = solver.Acidogenesis_LIN(alpha, D, S1, X1, XT_min.x, S1in, XT)
    
    k6 = CH4_min.x
    khyd = XT_min.x  
    
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
    fprint('LIN', 'kLa', [kLa, -kH], r2_CO2, err_CO2)
    fprint('LIN', 'k1', k1, r2_S1, err_S1)
    
    # Visualization
    Y_model_CH4 = np.empty(len(D))
    Y_model_CO2 = np.empty(len(D))
    Y_data_CO2  = np.empty(len(D))
    Y_model_XT  = np.empty(len(D))
    Y_model_S1  = np.empty(len(D))
    Y_data_S1   = np.empty(len(D))
    Y_data_XT   = np.empty(len(D))
    X_data_CO2  = np.empty(len(D))
    X_data_S1   = np.empty(len(D))
    
    for i in range(0, len(D)):
        # Methane Flowrate
        Y_model_CH4[i] = k6 * alpha * D[i] * X2[i]
        # Carbon Dioxide Flowrate
        pKb = -np.log10(Kb)
        X_data_CO2[i] = C[i] / (1 + 10**(pH[i] - pKb)) / PC[i] 
        Y_model_CO2[i] = kLa * X_data_CO2[i] + kH
        Y_data_CO2[i]  = qC[i]/PC[i]
        # Hydrolysis
        Y_data_XT[i] = D[i] * (XTin - XT[i])
        Y_model_XT[i] = khyd * XT[i]    
        # Acidogenesis
        Y_data_S1[i] = D[i]*(S1in - S1[i]) + khyd * XT[i] 
        X_data_S1[i] = alpha * D[i] * X1[i]
        Y_model_S1[i] = k1 * X_data_S1[i]

    regplot(X2, [qM, Y_model_CH4], 'X2 [mol/m3]', 'qM [mol/m3/d]', 'k6 regression')
    regplot(X_data_CO2, [Y_data_CO2, Y_model_CO2], 'C/(1+10**(pH-pKb))-KH*PC [mol/m3]', 'qC [mol/m3/d]', 'kLa regression')
    regplot(XT, [Y_data_XT, Y_model_XT], 'XT [kg/m3]', 'D*(XTin-XT) [kg/m3/d]', 'khyd regression')
    regplot(X_data_S1, [Y_data_S1, Y_model_S1], 'alpha*D*X1 [kg/m3/d]', 'D*(S1in-S1)+khyd*XT [kg/m3/d]', 'k1 regression', showbool=True)

    return 0;

if __name__ == "__main__":
    main()
    
