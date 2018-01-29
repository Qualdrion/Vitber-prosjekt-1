import numpy as np
from math import *

def fredholm_rhs (xc, F):
    '''Set up the RHS of the system
    INPUT :
        xc : defines collocation points
        F : function defining the geological survey measurements
    OUTPUT:
        vector defining the RHS of the system'''
    Nc = xc.shape [0]
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
    Nc = xc.shape[0]
    Ns = xs.shape[0]
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
        a.append(0.5 + 0.5*cos())
    return
