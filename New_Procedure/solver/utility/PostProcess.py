# utility function for printing and plotting
import matplotlib.pyplot as plt
import numpy as np
from typing import List
from sklearn.metrics import r2_score

def fprint(type: str, param: str, sol, score=0, flag='False', itr=0, funcall=0, intercept=0):
    if type == 'RMSE':
        str = f"== Results {param} w {type} =="
        print(str)
        print(f"Success: {flag}")
        print(f"Solution: {sol}")
        print(f"RMSE: {funcall}")
        print(f"Iterations: {itr}")
        print("="*len(str), '\n')
        
    elif type == 'LIN':
        str = f"== Results {param} w {type} =="
        print(str)
        print(f"Solution: {sol}")
        print(f"q: {intercept}")
        print(f"SCORE: {score}")
        print("="*len(str), '\n')
    
def modelplot(x, y: List, xlabel, ylabel, title, legend: List, showbool=False):
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
    
def regplot(x, y, xlabel, ylabel, title, showbool=False, gridbool=True, subplot=False):
    plt.figure(figsize=(10,5))
    
    if subplot == False:
        plt.plot(x, y[0], 'o', label = 'Experimental')
        plt.plot(x, y[1], '-', label = 'Model', linewidth=2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(gridbool)
        plt.legend()
        plt.title(title)
        
    elif subplot == True:
        plt.subplot(1,2,1)
        plt.plot(x[0], y[0], 'o', label = 'Experimental')
        plt.plot(x[0], y[1], '-', label = 'Model', linewidth=2)
        plt.xlabel(xlabel[0])
        plt.ylabel(ylabel[0])
        plt.grid(gridbool)
        plt.legend()
        
        plt.subplot(1,2,2)
        plt.plot(x[1], y[2], 'o', label = 'Experimental')
        plt.plot(x[1], y[3], '-', label = 'Model', linewidth=2)
        plt.xlabel(xlabel[1])
        plt.ylabel(ylabel[1])
        plt.grid(gridbool)
        plt.legend()
        plt.suptitle(title)
    
    if showbool == True:
        plt.show()
    else:
        pass
    
def multiplot(z:List, x: List, xlabel: List, zlabel, title, showbool=False, gridbool=True):
    fig = plt.figure(figsize=(10,5))
    ax  = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_title(title)
    ax.set_xlabel(xlabel[0])
    ax.set_ylabel(xlabel[1])
    ax.set_zlabel(zlabel)
    ax.plot(x[0], x[1], z[1], 'o', label = 'Model')
    ax.plot(x[0], x[1], z[0], '-', label = 'Experimental', linewidth=2)
    
    if showbool == True:
        plt.show()
    else:
        pass