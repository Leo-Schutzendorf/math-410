import numpy as np
import math
import matplotlib.pyplot as plt
import scipy

def fixedpoint(N,g,p_0,TOL):
    iterListx = []
    iterListy = []
    i=1
    p_prev=p_0#new variable not imported into function
    while i<=N:
        p=g(p_prev)
        if abs(p-p_prev)<TOL:
            return p,iterListx,iterListy

        print("$p_"+str(i),"=",str(p)+"$\n")
        i+=1
        p_prev=p
        iterListx.append(p)
        iterListy.append(g(p))
    print("The method failed after",str(i),"iterations")
    return p,iterListx,iterListy

N=144
g=lambda x:np.exp(x)-x-1
p_0=1
TOL = 0.01

answer,xlist,ylist = fixedpoint(N,g,p_0,TOL)

print(answer)

x=np.linspace(-2, 2, 1000)
y=g(x)
plt.plot(np.array(xlist),np.array(ylist))
plt.show()