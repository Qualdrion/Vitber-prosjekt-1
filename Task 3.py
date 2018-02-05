import numpy as np
from math import *
from matplotlib import pyplot as plt
import pickle
import scipy.special as ss
import sympy as sy
from Test_Example import analytical_solution

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
        b[i] = F(xc[i], 0.025)
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
                A[i][j] += K(xc[i], xq[k]) * w[k] * Lagrange_Basis(j, xq[k], xs, Ns)
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
            w.append(0.5/n)
        else:
            w.append(1/n)
    return xq, w

def Legendre(n):
    x1, w1 = np.polynomial.legendre.leggauss(n)
    xq = [0.5*(x+1) for x in x1]
    w = [0.5*x for x in w1]
    return xq, w

def Lagrange_Basis (j, xq, xs, ran):
    L = 1
    for i in range(ran):
        if j != i:
            L *= (xq-xs[i])/(xs[j]-xs[i])
    return L

def Density(a):
    p = []
    for i in range(len(a)):
        p.append(sin(3*pi*a[i])*exp(-2*a[i]))
    return p

def Gen_Error(p, xc, xs, K, F, method):
    x = []
    y = []
    b = fredholm_rhs(xc, F)
    for i in range(1, 9):
        print(i)
        xq, w = method(2**i)
        x.append(2**i)
        A = fredholm_lhs(xc, xs, xq, w, K)
        Ap = np.dot(A, p).tolist()
        r = []
        for j in range(len(Ap)):
            r.append(abs(Ap[j]-b[j]))
        y.append(max(r))
    return x, y

def Gen_Error_p(start, end, K, F, method):
    x = []
    y = []
    for i in range(start, end+1):
        print(i)
        xc = Chebyshev(i)
        b = fredholm_rhs(xc, F)
        xq, w = method(i**2)
        A = fredholm_lhs(xc, xc, xq, w, K)
        p = np.linalg.solve(A, b)
        p2 = Density(xc)
        r = []
        for j in range(len(p)):
            r.append(abs(p[j]-p2[j]))
        x.append(i)
        y.append(max(r))
    return x, y

def Gen_Perturbed(n, F):
    x = []
    y = []
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    plt.plot(xc, b)
    plt.plot(xc, b2)
    plt.show()

def Gen_plot_perturbed(n, F, K, method):
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    xq, w = method(10*n)
    A = fredholm_lhs(xc, xc, xq, w, K)
    p1 = np.linalg.solve(A, b)
    p2 = np.linalg.solve(A, b2)
    p3 = Density(xc)
    plt.plot(xc, p1)
    plt.plot(xc, p2)
    plt.plot(xc, p3)
    plt.show()

def Tikhonov(n, F, K, method):
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    xq, w = method(10*n)
    A = fredholm_lhs(xc, xc, xq, w, K)
    p1 = Density(xc)
    e = []
    l = []
    for i in range(-14, 2):
        p = []
        print(i)
        lhs = np.dot(A.T, A) + np.dot(10**i, np.identity(n))
        rhs = np.dot(A.T, b)
        p.append(np.linalg.solve(lhs, rhs))
        r = []
        for j in range(len(p)):
            r.append(abs(p[j]-p1[j]))
        e.append(max(r[0].tolist()))
        l.append(10**i)
    plt.plot(l, e)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

def Plot_func(x, y):
    for yi in y:
        plt.plot(x, yi)
    plt.yscale('log')
    plt.show()


'''a = Chebyshev(40)
p = Density(a)
K = lambda x, y: 0.025 * (0.025**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
x, y1 = Gen_Error(p, a, a, K, F, Legendre)
x, y2 = Gen_Error(p, a, a, K, F, Trapezoid)
y = [y1, y2]
Plot_func(x, y)'''


'''K1 = lambda x, y: 0.025 * (0.025**2 + (y-x)**2)**(-3/2)
K2 = lambda x, y: 0.25 * (0.25**2 + (y-x)**2)**(-3/2)
K3 = lambda x, y: 2.5 * (2.5**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
x, y1 = Gen_Error_p(5, 30, K1, F, Legendre)
x, y2 = Gen_Error_p(5, 30, K2, F, Legendre)
x, y3 = Gen_Error_p(5, 30, K3, F, Legendre)
y = [y1, y2, y3]
Plot_func(x, y)'''

'''F = pickle.load( open( "F.pkl", "rb" ) )
Gen_Perturbed(30, F)'''

K1 = lambda x, y: 0.025 * (0.025**2 + (y-x)**2)**(-3/2)
K2 = lambda x, y: 0.25 * (0.25**2 + (y-x)**2)**(-3/2)
K3 = lambda x, y: 2.5 * (2.5**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
#Gen_plot_perturbed(30, F, K3, Legendre)
Tikhonov(30, F, K3, Legendre)