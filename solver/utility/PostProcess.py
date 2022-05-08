# utility function for printing and plotting
import matplotlib.pyplot as plt
import numpy as np
from typing import List
from sklearn.metrics import r2_score

def fprint(type: str, param: str, sol, r2=0, err=0, flag='False', itr=0, funcall=0):
    if type is 'RMSE':
        str = f"== Results {param} w {type} =="
        print(str)
        print(f"Success: {flag}")
        print(f"Solution: {sol}")
        print(f"RMSE: {funcall}")
        print(f"Iterations: {itr}")
        print("="*len(str), '\n')
        
    elif type is 'LIN':
        str = f"== Results {param} w {type} =="
        print(str)
        print(f"Solution: {sol}")
        print(f"R2: {r2}")
        print(f"ERR: {err}")
        print("="*len(str), '\n')
    
def regplot(x, y: List, xlabel, ylabel, title, legend: List, showbool=False):
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