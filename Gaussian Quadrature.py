import numpy as np
import scipy.integrate
import math

def GaussQuad(n, f):
    # Retrieve the Gauss-Legendre roots (nodes) and corresponding weights
    roots, weights = np.polynomial.legendre.leggauss(n)
    
    # Vectorize the scalar math function to handle NumPy array inputs safely
    f_vec = np.vectorize(f)
    
    # Calculate the Gauss terms and return their sum
    return np.sum(weights * f_vec(roots))

# Define order and target integration function
a=1
b=2
scale = (b-a)/2
n = 5
f = lambda x: math.cos((1/8)*math.pi*x + (1/8)*math.pi)**2

# Execute and print result
print(f"Approximated Integral: {GaussQuad(n, f)}")
