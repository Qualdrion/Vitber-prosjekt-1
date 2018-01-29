import numpy as np
from math import *
from matplotlib import pyplot as plt

def fredholm_rhs (xc, F):
    '''Set up the RHS of the system
    INPUT :
        xc : defines collocation points
        F : function defining the geological survey measurements
    OUTPUT:
        vector defining the RHS of the system'''
    Nc = len(xc)
    b = np.zeros(Nc)
    for i in range(Nc):
        b[i] = F(xc[i])
    return b

def fredholm_lhs (xc, xs, xq, w, K):
    '''
    Set up the LHS of the system
    INPUT:
    xc: collocation points
    xs: source points
    xq, w: numerical quadrature
    K: integral kernel
    OUTPUT:
    matrix defining the LHS of the system'''
    Nc = len(xc)
    Ns = len(xs)
    A = np.zeros((Nc, Ns))
    #FIXIT : implement the function!
    for i in range(Nc):
        for j in range(Ns):
            for k in range(len(xq)):
                L = 1
                for m in range(Ns):
                    if m != j:
                        L = L * (xq[k]-xs[m])/(xs[j]-xs[m])
                A[i][j] += K(xc[i], xq[k]) * L * w[k]
    return A

def Chebyshev(n):
    a = []
    for i in range(1, n+1):
        a.append(0.5 + 0.5*cos(pi*(2*i-1)/(2*n)))
    return a

def Trapezoid(n):
    xq = []
    w = []
    for i in range(n+1):
        xq.append(i/n)
        if i == 0 or i == n:
            w.append(0.5)
        else:
            w.append(1)
    return xq, w

def Density(a):
    p = []
    for i in range(len(a)):
        p.append(sin(3*pi*a[i])*exp(-2*a[i]))
    return p

def Gen_Error(n, p, xc, xs, K, F):
    y = []
    for i in range(1, 8):
        xq, w = Trapezoid(2**i)
        A = fredholm_lhs(xc, xs, xq, w, K)
        b = fredholm_rhs(xc, F)
        Ap = np.dot(A, p).tolist()
        r = []
        for j in range(len(Ap)):
            r.append(Ap[j]-b[j])
        y.append(max(r))
    return y


def Plot_func(y):
    x = np.arange(1, len(y)+1)
    plt.plot(x, y)
    plt.show()

a = Chebyshev(40)
p = Density(a)
K = lambda x, y: 0.025 * (0.025**2 + (y-x)**2)**(-3/2)
F = lambda x: (1-x)/(0.025*(0.025**2+(x-1)**2)**(1/2)) - (-x)/(0.025*(0.025**2+(x)**2)**(1/2))
y = Gen_Error(100, p, a, a, K, F)
Plot_func(y)