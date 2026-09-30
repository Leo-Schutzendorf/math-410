import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

def bisecshit(a,b,TOL,mathFunc,N):
    #a is left bound
    #b is right bound
    #TOL is epsilon
    #mathFunc is the mathematical function
    i=1 # counter

    if (mathFunc(a)<0 and mathFunc(b)<0) or (mathFunc(a)>0 and mathFunc(b)>0):
        raise Exception("f(a) and f(b) must have different signs")
    
    FA = mathFunc(a)
    p=0.8
    p_old=p+67/TOL#p_{n-1}

    print("Press 1 to do N steps in the while loop, like in Algorithm 2.1")
    print("Press 2 to end the while loop once the relative error goes under the tolerance value")
    method = input("Which method?")
    if method[0] == "1" or method == "One" or method == "one":
        while i<=N:
            if i>1:
                p=a+(b-a)/2
            p_old=p#p_{n-1}

            FP = float(mathFunc(p))
            if FP==0 or (b-a)/2<TOL:
                return(p)
            if FA*FP>0:
                a=p
                FA=FP
            else:
                b=p
            print("p_{"+str(i)+"} =",str(p) +",","New a =",str(a) +",","New b =",str(b))
            print("\n")
            i+=1
        raise Exception("Method failed after N iterations")
    else:
        while i<=N:
            p_old=p#p_{n-1}

            FP = float(mathFunc(p))
            if FP==0 or abs((p-p_old)/p)>=TOL:
                return(p)
            if FA*FP>0:
                a=p
                FA=FP
            else:
                b=p
            print("p_{"+str(i)+"} =",str(p) +",","New a =",str(a) +",","New b =",str(b))
            print("\n")
            p=a+(b-a)/2
            
            print("Current progress")
            print("On iteration #"+str(i))
            i+=1
        return p


a=0.1
b=0.81
TOL = 0.0001
mathFunc = lambda x:((2-3*x*x)**0.25)-x
itmax = 4 #aka N

if a>b:
    a,b = b,a
answer = bisecshit(a,b,TOL,mathFunc,itmax)
print("Approximated root: x =",str(answer),"and f(x) =",str(mathFunc(answer)))
print("Upper bound on error:",str((b-a)/2**itmax))

x=np.linspace(0.75, 2.25, 1000)
y=mathFunc(x)
plt.plot(x,y)
plt.show()