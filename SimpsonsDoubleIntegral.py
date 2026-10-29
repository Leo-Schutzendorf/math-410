import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

def simpsonsDouble(a, b, m, n, c, d, f):
    h = (b - a) / n
    J1 = 0.0   # end terms   (i = 0 and i = n)
    J2 = 0.0   # even terms
    J3 = 0.0   # odd  terms

    for i in range(n + 1):
        x = a + i * h
        HX = (d(x) - c(x)) / m          # <-- fixed: d - c

        # endpoints only
        K1 = f(x, c(x)) + f(x, d(x))
        K2 = 0.0   # even interior
        K3 = 0.0   # odd  interior

        for j in range(1, m):           # interior points only
            y = c(x) + j * HX
            Q = f(x, y)
            if j % 2 == 0:
                K2 += Q
            else:
                K3 += Q

        # Composite Simpson in y
        L = (K1 + 2 * K2 + 4 * K3) * HX / 3

        if i == 0 or i == n:
            J1 += L
        elif i % 2 == 0:                # <-- fixed: test i, not J2
            J2 += L
        else:
            J3 += L

    return h * (J1 + 2 * J2 + 4 * J3) / 3

''' test simpsonsDouble
a, b = 0,1
c = lambda x: (2+x)/3
d = lambda x: x+x
n, m = 2,2
f = lambda x, y: (x**2)+(y**3)
print(simpsonsDouble(a, b, m, n, c, d, f))
'''

def GaussianTriple(a, b, m, n, p, c, d, alpha, beta, f):
    # Generate Gaussian roots and coefficients for each order
    rx, cx = np.polynomial.legendre.leggauss(n)
    ry, cy = np.polynomial.legendre.leggauss(m)
    rz, cz = np.polynomial.legendre.leggauss(p)

    J = 0.0

    # Integrate in x
    for i in range(n):
        x = (b - a) * rx[i] / 2 + (b + a) / 2
        wx = cx[i] * (b - a) / 2

        # Integrate in y
        hy = (d(x) - c(x)) / 2
        my = (d(x) + c(x)) / 2
        K = 0.0

        for j in range(m):
            y = hy * ry[j] + my
            wy = cy[j] * hy

            # Integrate in z
            hz = (beta(x, y) - alpha(x, y)) / 2
            mz = (beta(x, y) + alpha(x, y)) / 2
            L = 0.0

            for k in range(p):
                z = hz * rz[k] + mz
                wz = cz[k] * hz

                L += wz * f(x, y, z)

            K += wy * L

        J += wx * K

    return J


# Test
a, b = 0, 1
c = lambda x: 1
d = lambda x: x
alpha = lambda x, y: 0
beta = lambda x, y: y

n, m, p = 2, 2, 2
f = lambda x, y, z: y*y*z

print(GaussianTriple(a, b, m, n, p, c, d, alpha, beta, f))
