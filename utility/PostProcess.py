# utility function for printing and plotting
import matplotlib.pyplot as plt
import numpy as np
from typing import List
from sklearn.metrics import r2_score

def fprint(type, sol, flag, itr, funcall, r2: List):
    str = f"== Results {type} =="
    print(str)
    print(f"Success: {flag}")
    print(f"Solution: {sol}")
    print(f"RSE: {funcall}")
    print(f"Iterations: {itr}")
    print(f"R2 score: {r2_score(r2[0], r2[1])}")
    print("="*len(str), '\n')
    
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