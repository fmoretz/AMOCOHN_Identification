# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve
from scipy.optimize import Bounds
from typing import List
from sklearn.metrics import r2_score
plt.close('all')

# import data and define the variables
global C_in, q_C, X_1, X_2, C, q_C, D, kd1, kd2
from dataimport import *
from Identification import *
kd1 = kd[0]; kd2 = kd[1]
C_in  = float(T2["Cin"])
D     = 1/T1["HRT"]
X_1   = T1["X1"]
X_2   = T1["X2"]
C     = T1["C"]
q_C   = T1["q_C"]        # Experimental points for the regression and comparison

# create main function
def main():
    
    # define the bounds
    bd = Bounds([0,0],[1000,1000])
    
    # define the initial values
    gs = [1, 10]
    
    # solve minimization problem
    minsol = minimize(rmse, gs, args = (D, X_1, X_2, C, q_C, C_in, kd1, kd2, 1), method = 'SLSQP', bounds = bd, tol = 1e-10)
    qCmodel = qCeval(D, X_1, X_2, C, C_in, kd1, kd2, alfa, minsol.x[0], minsol.x[1])
    
    # print out results
    fprint(minsol.x, minsol.success, minsol.nit, minsol.fun, [q_C, qCmodel])

    # plot results
    regplot(1/D, [q_C, qCmodel], 'time [days]', 'q_C [g/s]', 'Carbon dioxide production regression', ['Experimental data', 'Model'], False)
    
    # plot function surface
    k4x = np.linspace(0, 120, 20)   
    k5x = np.linspace(0, 120, 20)
    
    rmsez = np.empty([len(k4x), len(k5x)])
    for i in range(0,len(k4x)):
        for j in range(0,len(k5x)):
            rmsez[i][j] = rmse([k4x[i], k5x[j]], D, X_1, X_2, C, q_C, C_in, kd1, kd2, 1)
    
    surfplot(k4x, k5x, rmsez, 'k4', 'k5', 'rmse', 'rmse over k4 and k5 values', True)

    return 0;


# define the function for RSE
def rmse(x, D, X1, X2, C, qC, Cin, kd1, kd2, alfa) -> float:
    qClist = []
    for i in range(0, len(D)):
         qClist.append(qCeval(D[i], X1[i], X2[i], C[i], Cin, kd1, kd2, alfa, x[0], x[1]))
         
    return np.sqrt(np.mean(qC-qClist)**2)
    
# define the function for the qC
def qCeval(D, X1, X2, C, Cin, kd1, kd2, alfa, k4, k5) -> float:
    return D*(Cin-C) + k4*(alfa*D+kd1)*X1 + k5*(alfa*D+kd2)*X2

# utility function for printing and plotting
def fprint(sol, flag, itr, funcall, r2: List):
    print(f"\nSuccess: {flag}\n")
    print("==Results==")
    print(f"Solution: {sol}")
    print(f"RSE: {funcall}")
    print(f"Iterations: {itr}")
    print(f"R2 score: {r2_score(r2[0], r2[1])}")
    print("============")
    
def regplot(x, y: List, xlabel, ylabel, title, legend: List, showbool):
    plt.figure(figsize=(10,5))
    plt.plot(x, y[0], 'o', label = legend[0])
    plt.plot(x, y[1], '-', label = legend[1], linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.legend()
    plt.title(title)
    
    if showbool == True:
        plt.show()
    else:
        pass

def surfplot(x, y, z, xlabel, ylabel, zlabel, title, showbool):
    import mpl_toolkits.mplot3d.axes3d as axes3d
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    
    if showbool == True:
        plt.show()
    else:
        pass
    
    

if __name__ == "__main__":
    main()
    
