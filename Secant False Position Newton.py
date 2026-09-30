import math
import numpy as np
import matplotlib.pyplot as plt
import scipy

def find_root(p_0,p_1,N,TOL,f):
    whichone = input("Which method?  Press S for Secant, F for False Position, or N for Newton's.")
    if whichone[0] == "S" or whichone[0] == "s":#Secant
        i=2
        q_0 = f(p_0)
        q_1 = f(p_1)
        while i<=N:
            p = p_1-(q_1*(p_1-p_0)/(q_1-q_0))
            if abs(p-p_1)<TOL:
                return p

            i+=1

            p_0=p_1
            q_0=q_1
            p_1=p
            q_1=f(p)

            if i-2==1:
                print(f"Current estimate is {p} after {i-2} iteration.")
            else:
                print(f"Current estimate is {p} after {i-2} iterations.")
        print(f"The method failed after {N} iterations.")
        return(p)
    
    if whichone[0] == "F" or whichone[0] == "f":#False position
        i=2
        q_0 = f(p_0)
        q_1 = f(p_1)
        while i<=N:
            p = p_1-(q_1*(p_1-p_0)/(q_1-q_0))
            if abs(p-p_1)<TOL:
                return p

            i+=1
            q=f(p)

            if q*q_1<0:
                p_0=p_1
                q_0=q_1
            
            p_1=p
            q_1=q

            if i-2==1:
                print(f"Current estimate is {p} after {i-2} iteration.")
            else:
                print(f"Current estimate is {p} after {i-2} iterations.")
            
        print(f"The method failed after {N} iterations.")
        return(p)
    if whichone[0] == "N" or whichone[0] == "n":  # Newton's
        df = lambda x: float(scipy.differentiate.derivative(f, x).df)
        i = 1

        while i <= N:
            slope = df(p_0)
            if slope == 0:
                print("Derivative is already 0.")
                return p_0

            p = p_0 - f(p_0) / slope

            if abs(p - p_0) < TOL:
                return p

            if i==1:
                print(f"Current estimate is {p} after {i} iteration.")
            else:
                print(f"Current estimate is {p} after {i} iterations.")
            i += 1
            p_0 = p

        print(f"The method failed after {N} iterations.")
        return p
f = lambda x:np.exp(6*x) + 3*(math.log(2)**2)*np.exp(2*x) - math.log(8)*np.exp(4*x) - (math.log(2)**3)


N=73
p_0=0
p_1=2
TOL = 0.0002
print("Answer = "+str(find_root(p_0,p_1,N,TOL,f)))