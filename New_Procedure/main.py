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
    
    # Single Parameter Regression
    [k6, q_ME, score_ME] = solver.Methane_LIN(alpha, D, X2, qM)
    [khyd, q_XT, score_XT]	= solver.Hydrolysis_LIN(D, XTin, XT)
    [k1, q_AC, score_AC] = solver.Acidogenesis_LIN(alpha, D, S1, X1, khyd, S1in, XT)
    [kLa, q_CO2, score_CO2] = solver.Carbon_Dioxide_LIN(CO2, PC, KH, qC)
    [[a, b], q_1, score_1] = solver.Substrate1_LIN(S1, alpha, D)
    [[c, d, e, f], q_2, score_2] = solver.Substrate2_LIN(S2, alpha, D)
    [[k6_k3, k6k2_k3k1], q_3, score_3] = solver.Methanogenesis_RATIO(D, qM, S2in, S2, S1in, S1, khyd, XT)
    [[k4_k1, k5_k6], q_4, score_4] = solver.Carbogenesis_RATIO(D, qC, qM, Cin, C, S1in, S1, khyd, XT)
    
    # Parameter retrieving
    mu1m = 1/a
    Ks1  = q_1/0.11
    mu2m = 1/c
    Ks2  = q_2/0.11
    KI   = 1/(e*mu2m)
    
    k3   = k6 / k6_k3
    k2   = k6k2_k3k1 / k6_k3 * k1
    k4   = k4_k1 * k1
    k5   = k5_k6 * k6 
    
    # Model evaluation
    qM_Model = np.empty(len(D))
    XT_Model = np.empty(len(D))
    
    for i in range(0, len(D)):
        qM_Model[i] = rm.qMeval(k6, alpha, D[i], X2[i])
        XT_Model[i] = rm.XTeval(khyd, D[i], XTin)

    # print out results
    fprint('LIN', 'k6', k6, score=score_ME, intercept=q_ME)
    fprint('LIN', 'khyd', khyd, score=score_XT, intercept=q_XT)
    fprint('LIN', 'k1', k1, score=score_AC, intercept=q_AC)
    fprint('LIN', 'kLa', kLa, score=score_CO2)
    fprint('LIN', 'mu1m + Ks1', [mu1m, Ks1], score=score_1, intercept=q_1)
    fprint('LIN', 'mu2m + Ks2 + KI', [mu2m, Ks2, KI], score=score_2, intercept=q_2)
    fprint('LIN', 'k3 + k2', [k3, k2], score=score_3, intercept=q_3)
    fprint('LIN', 'k5 + k4', [k5, k4], score=score_4, intercept=q_4)

    
    
    # Visualization
    [X_CH4, Y_data_CH4, Y_model_CH4] = rm.CH4model(alpha, D, X2, qM, k6)
    [X_XT, Y_data_XT, Y_model_XT] = rm.XTmodel(D, XTin, XT, khyd)
    [X_CO2, Y_data_CO2, Y_model_CO2] = rm.CO2model(qC, CO2, KH, PC, kLa)
    [X_AC, Y_data_AC, Y_model_AC] = rm.ACmodel(alpha, D, S1in, S1, XT, X1, k1, khyd)
    [X_S1_1, X_S1_2, Y_data_S1, Y_model_S1] = rm.S1model(alpha, D, mu1m, Ks1, S1)
    [X_ME_1, X_ME_2, Y_data_ME, Y_model_ME] = rm.MEmodel(D, qM, S2in, S2, S1in, S1, khyd, XT, k6, k3, k1, k2)
    [X_IC_1, X_IC_2, Y_data_IC, Y_model_IC] = rm.ICmodel(D, qC, qM, Cin, C, S1in, S1, khyd, XT, k4, k1, k5, k6)
    
    regplot(X_CH4,  [Y_data_CH4, Y_model_CH4], 'D [1/d]', 'qM/X2 [1/d]', 'k6 regression - Methane Flowrate')
    regplot(X_CO2,  [Y_data_CO2, Y_model_CO2], 'CO2-KH*PC [mol/m3]', 'qC [mol/m3/d]', 'kLa regression - Carbon Dioxide Flowrate')
    regplot(X_XT,   [Y_data_XT,  Y_model_XT],  'XT [kg/m3]', 'D*(XTin-XT) [kg/m3/d]', 'khyd regression - Hydrolysis')
    regplot(X_AC,   [Y_data_AC,  Y_model_AC],  'alpha*D*X1 [kg/m3/d]', 'D*(S1in-S1)+khyd*XT [kg/m3/d]', 'k1 regression - Acidogenesis')
    regplot(
            x=[X_S1_1, X_S1_2], 
            y=[Y_data_S1, Y_model_S1, Y_data_S1, Y_model_S1], 
            xlabel=['alpha*D*S1/0.9 [kg/m3/d]', 'alpha*D/0.9 [d-1]'],
            ylabel=['S1 [kg/m3]', 'S1 [kg/m3]'],
            title='mu1m & Ks1regression - Substrate 1',
            subplot=True
            )
    regplot(
            x=[X_ME_1, X_ME_2], 
            y=[Y_data_ME, Y_model_ME, Y_data_ME, Y_model_ME], 
            xlabel=['alpha*D*X2 [kg/m3/d]', '-alpha*D*X1 [kg/m3/d]'],
            ylabel=['D*(S2in-S2) [kg/m3/d]', 'D*(S2in-S2) [kg/m3/d]'],
            title='k2 & k3 regression - Methanogenesis',
            subplot=True
            )
    regplot(
            x=[X_IC_1, X_IC_2], 
            y=[Y_data_IC, Y_model_IC, Y_data_IC, Y_model_IC], 
            xlabel=['alpha*D*X2 [kg/m3/d]', '-alpha*D*X1 [kg/m3/d]'],
            ylabel=['qC-D*(Cin-C) [kg/m3/d]', 'qC-D*(Cin-C) [kg/m3/d]'],
            title='k4 & k5 regression - Inorganic Carbon',
            showbool=True,
            subplot=True
            )
    
    multiplot([Y_data_S1, Y_model_S1], [X_S1_1, X_S1_2], ['alpha*D*S1/0.9 [kg/m3/d]','alpha*D/0.9 [d-1]'], 'S1 [kg/m3]', 'mu1m and Ks1 regression - Substrate 1')

    return 0;

if __name__ == "__main__":
    main()
    
