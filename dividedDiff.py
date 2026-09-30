import math
import numpy as np
import matplotlib.pyplot as plt
import scipy

def divided_diff(xs, ys):
    
    n = len(ys)
    fTable = np.zeros([n, n])
    fTable[:, 0] = ys
    
    for j in range(1, n):
        for i in range(n - j):
            fTable[i][j] = (fTable[i + 1][j - 1] - fTable[i][j - 1]) / (xs[i + j] - xs[i])

    print("F-table:",fTable)
    return fTable[0, :]

def newton_poly(coef, x_data, x):
    n = len(x_data) - 1
    p = coef[n]
    for k in range(1, n + 1):
        p = coef[n - k] + (x - x_data[n - k]) * p
    return p

x_data = np.array([-1/4,1/4])
y_data = np.array([1.33203, 0.800781])

coefficients = divided_diff(x_data, y_data)
print("Newton Coefficients:", coefficients)

x_val = 0.43
y_val = newton_poly(coefficients, x_data, x_val)
print(f"Interpolated value at x = {x_val} is y = {y_val:.4f}")