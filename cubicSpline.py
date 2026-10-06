import numpy as np

def cubicSpline(x_s, y_s,FPO=None,FPN=None):
    whichOne = input("Natural or clamped?")
    if whichOne == "n" or whichOne == "N" or whichOne == "natural" or whichOne == "Natural" or FPO ==None or FPN == None:
        # n is the number of intervals, which is len(x_s) - 1
        n = len(x_s) - 1
        
        # Step 1: Calculate interval widths (h_i)
        h_s = [x_s[i+1] - x_s[i] for i in range(n)]
        
        # Step 2: Compute alpha values
        alpha = [0] * n
        for i in range(1, n):
            alpha[i] = (3 / h_s[i]) * (y_s[i+1] - y_s[i]) - (3 / h_s[i-1]) * (y_s[i] - y_s[i-1])
            
        # Steps 3 & 4: Forward substitution to solve the tridiagonal system (l, mu, z)
        l_s = [1] * (n + 1)
        mu_s = [0] * (n + 1)
        z_s = [0] * (n + 1)
        
        for i in range(1, n):
            l_s[i] = 2 * (x_s[i+1] - x_s[i-1]) - h_s[i-1] * mu_s[i-1]
            mu_s[i] = h_s[i] / l_s[i]
            z_s[i] = (alpha[i] - h_s[i-1] * z_s[i-1]) / l_s[i]
            
        # Step 5: Boundary conditions for natural cubic spline
        l_s[n] = 1
        z_s[n] = 0
        c_s = [0] * (n + 1)
        b_s = [0] * n
        d_s = [0] * n
        
        # Step 6: Backward substitution to find coefficients b, c, d
        for j in range(n - 1, -1, -1):
            c_s[j] = z_s[j] - mu_s[j] * c_s[j+1]
            b_s[j] = (y_s[j+1] - y_s[j]) / h_s[j] - h_s[j] * (c_s[j+1] + 2 * c_s[j]) / 3
            d_s[j] = (c_s[j+1] - c_s[j]) / (3 * h_s[j])
            
        # Return coefficients (a, b, c, d) for each of the n pieces
        return [(y_s[j], b_s[j], c_s[j], d_s[j]) for j in range(n)]
    else:
        n = len(x_s) - 1

        # Step 1: Calculate interval widths (h_i)
        h_s = [x_s[idx+1] - x_s[idx] for idx in range(n)]
        
        # Step 2 & 3: Compute alpha values for interior points and boundaries
        alpha_s = [0.0] * (n + 1)
        
        # Safely extract first and last elements using explicit variables
        h0 = h_s[0]
        hn_minus_1 = h_s[n-1]
        y0 = y_s[0]
        y1 = y_s[1]
        yn = y_s[n]
        yn_minus_1 = y_s[n-1]
        
        alpha_s[0] = (3.0 * (y1 - y0) / h0) - 3.0 * FPO
        alpha_s[n] = 3.0 * FPN - (3.0 * (yn - yn_minus_1) / hn_minus_1)
        
        for i in range(1, n):
            alpha_s[i] = (3.0 / h_s[i]) * (y_s[i+1] - y_s[i]) - (3.0 / h_s[i-1]) * (y_s[i] - y_s[i-1])
            
        # Step 4: Forward substitution initialization
        l_s = [0.0] * (n + 1)
        mu_s = [0.0] * (n + 1)
        z_s = [0.0] * (n + 1)
        
        l_s[0] = 2.0 * h0
        mu_s[0] = 0.5
        z_s[0] = alpha_s[0] / l_s[0]
        
        # Step 5: Forward substitution for interior tridiagonal elements
        for i in range(1, n):
            l_s[i] = 2.0 * (x_s[i+1] - x_s[i-1]) - h_s[i-1] * mu_s[i-1]
            mu_s[i] = h_s[i] / l_s[i]
            z_s[i] = (alpha_s[i] - h_s[i-1] * z_s[i-1]) / l_s[i]
            
        # Step 6: Solve boundary conditions for the last node
        l_s[n] = hn_minus_1 * (2.0 - mu_s[n-1])
        z_s[n] = (alpha_s[n] - hn_minus_1 * z_s[n-1]) / l_s[n]
        
        # Initialize coefficients lists
        a_s = list(y_s)
        b_s = [0.0] * n
        c_s = [0.0] * (n + 1)
        d_s = [0.0] * n
        
        c_s[n] = z_s[n]
        
        # Step 7: Backward substitution to find coefficients b, c, d
        for j in range(n - 1, -1, -1):
            c_s[j] = z_s[j] - mu_s[j] * c_s[j+1]
            b_s[j] = (a_s[j+1] - a_s[j]) / h_s[j] - h_s[j] * (c_s[j+1] + 2.0 * c_s[j]) / 3.0
            d_s[j] = (c_s[j+1] - c_s[j]) / (3.0 * h_s[j])
            
        return [(a_s[j], b_s[j], c_s[j], d_s[j]) for j in range(n)]


x_s = [-0.25, 0.25]
y_s = [1.33203, 0.800781]
FPO = 2#f'(x_{0})
FPN = 5.43656#f'(x_{n})
print(cubicSpline(x_s, y_s,FPO,FPN))