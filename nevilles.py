import math
import numpy as np
import matplotlib.pyplot as plt
import scipy

def neville_interpolation(x, y, x_val):
    """
    Evaluates the interpolating polynomial at x_val using Neville's iterated method.
    
    Parameters:
    x (list): List of x-coordinates of the data points.
    y (list): List of y-coordinates of the data points.
    x_val (float): The value at which to evaluate the polynomial.
    
    Returns:
    float: The interpolated value at x_val.
    """
    n = len(x) - 1
    Q = [[0.0] * (n + 1) for _ in range(n + 1)]
    
    for i in range(n + 1):
        Q[i][0] = y[i]
        
    for j in range(1, n + 1):
        for i in range(n - j + 1):
            Q[i][j] = ((x_val - x[i + j]) * Q[i][j - 1] + (x[i] - x_val) * Q[i + 1][j - 1]) / (x[i] - x[i + j])
            
    return Q[0][n],np.array(Q)

# Example usage:
x_points = [-2,-1,0,1,2]
y_points = [3**x for x in x_points]
x_val = 0.5

result,allResults = neville_interpolation(x_points, y_points, x_val)
print(f"Interpolated value at x = {x_val} is {result:.5f}")
print(f"Here's the table:\n{allResults}")
