import numpy as np

def hermite(xs, ys, dfs):
    n = len(xs)
    size = 2 * n
    
    # Initialize the Q table with the correct size (2n x 2n)
    Q = np.zeros((size, size))
    zs = []
    
    # Populate z array and the first two columns of Q
    for i in range(n):
        zs.append(xs[i])
        zs.append(xs[i])
        Q[2*i][0] = ys[i]
        Q[(2*i)+1][0] = ys[i]
        Q[(2*i)+1][1] = dfs[i]
        if i != 0:
            Q[2*i][1] = (Q[2*i][0] - Q[(2*i)-1][0]) / (zs[2*i] - zs[(2*i)-1])
            
    # Compute the remaining divided differences
    for i in range(2, size):
        for j in range(2, i + 1):
            Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (zs[i] - zs[i-j])
            
    # Return exactly 2n coefficients
    return [float(Q[var][var]) for var in range(size)]

xs = [-1, -0.5, 0, 0.5]
ys = [0.86199480, 0.95802009, 1.0986123, 1.2943767]
dfs = [0.1553624, 0.23269654, 1/3, 0.45186776]

print(hermite(xs, ys, dfs))
