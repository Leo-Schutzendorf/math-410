import numpy as np
import math

def steffensen_root_finder(function, x0, tol, N=100):
    """
    Finds a root of f(x) = 0 using Steffensen's root-finding variant
    (Aitken's acceleration applied to f(x)).
    """
    p_0 = x0
    i=1

    f = lambda x:function(x)+x

    while i<N:
        p_1 = f(p_0)
        p_2 = f(p_1)

        p = p_0-(((p_1-p_0)**2)/(p_2-(2*p_1)+p_0))

        if abs(p-p_0)<tol:
            print("Converged after",str(i),"iterations.")
            return(p)

        i+=1
        p_0=p
        if i-1>=2:
            print("Estimate =",str(p),"after",str(i-1),"iterations") 
        else:
            print("Estimate =",str(p),"after",str(i-1),"iteration") 

# --- Your exact expression ---
def my_function(x):
    return (np.exp(6 * x) 
            + 3 * (math.log(2) ** 2) * np.exp(2 * x) 
            - math.log(8) * np.exp(4 * x) 
            - (math.log(2) ** 3))

# Because it contains exp(6*x), start close to 0 to avoid overflow
initial_guess = 0.0 

try:
    root = steffensen_root_finder(my_function, initial_guess,0.0002)
    print(f"Root found: {root:.10f}")
    print(f"f(root) evaluation: {my_function(root):.10e}")
except RuntimeError as e:
    print(e)
