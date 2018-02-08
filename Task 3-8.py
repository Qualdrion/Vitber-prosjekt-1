import numpy as np
from math import *
from matplotlib import pyplot as plt
import pickle
import scipy.special as ss
import sympy as sy
from Test_Example import analytical_solution
#The analytical solution of F comes from the code in the Test_Example.py file from the project website, and so
#The Test_Example file is imported such that we can use it to find F in the code.

#Method as specified in the assignment file. Has an extra input value d, as this is necessary for some of the problems.
def fredholm_rhs (xc, F, d):
    '''Set up the RHS of the system
    INPUT :
        xc : defines collocation points
        F : function defining the geological survey measurements
    OUTPUT:
        vector defining the RHS of the system'''
    Nc = len(xc)
    b = np.zeros(Nc)
    for i in range(Nc):
        b[i] = F(xc[i], d)
    return b


#Method as specified in the assignment file. Has an extra input value d, as this is necessary for some of the problems.
def fredholm_lhs (xc, xs, xq, w, K, d):
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
    #Triple for loop to generate the matrix A.
    for i in range(Nc):
        for j in range(Ns):
            for k in range(len(xq)):
                A[i][j] += K(xc[i], xq[k], d) * w[k] * Lagrange_Basis(j, xq[k], xs, Ns)
    return A

#Generates chebyshev's interpolation nodes with k as the input value n, in the interval [0, 1].
def Chebyshev(n):
    a = []
    for i in range(1, n+1):
        a.append(0.5 - 0.5*cos(pi*(2*i-1)/(2*n)))
    return a

#For the last problem it seemed like it might be relevant to have the option to generate Chebyshev's interpolation
#nodes in the general range a, b instead of the specific interval [0, 1], so this method is meant to do this.
def Chebyshev2(n, a, b):
    r = []
    for i in range(1, n+1):
        r.append(0.5*(a+b)+0.5*(a-b)*cos(pi*(2*i-1)/(2*n)))
    return r

#Returns nodes and weights for the trapezoid method in the interval [0, 1] for use in problem 3.
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

#Returns nodes and weights for the legendre gauss quadrature in the interval [0, 1] for use in problem 4-7.
def Legendre(n):
    x1, w1 = np.polynomial.legendre.leggauss(n)
    xq = [0.5*(x+1) for x in x1]
    w = [0.5*x for x in w1]
    return xq, w

#Returns nodes and weights for the legendre gauss quadrature in the generic interval [a, b] for use in problem 8.
def Legendre2(n, a, b):
    x1, w1 = np.polynomial.legendre.leggauss(n)
    xq = [((b-a)*x+b+a)/2 for x in x1]
    w = [0.5*(b-a)*x for x in w1]
    return xq, w

#Returns Lj(x) for use in creating the matrix A in most of the problems.
def Lagrange_Basis (j, xq, xs, ran):
    L = 1
    for i in range(ran):
        if j != i:
            L *= (xq-xs[i])/(xs[j]-xs[i])
    return L

#Returns the analytical density of p for all the x-values in the list a.
def Density(a):
    p = []
    for i in range(len(a)):
        p.append(sin(3*pi*a[i])*exp(-2*a[i]))
    return p

#For use in problem 3 and 4. Takes as input the list p with the densities at various points, xc and xs which are lists of collocation
#and source points, K and F which are functions, d, the depth, and method, which is either Legendre or Trapezoid depending on
#whether we want the graphs for problem 3 or 4.
def Gen_Error(p, xc, xs, K, F, method, d):
    x = []
    y = []
    b = fredholm_rhs(xc, F, d)
    #This for loop loops through Nq = 2^i for i=1, .., 8 and then appends 2^i, inf norm given Nq = 2^i to the lists x and y
    #And returns these 2 lists.
    for i in range(1, 20):
        print(i)
        xq, w = method(10*i)
        x.append(10*i)
        A = fredholm_lhs(xc, xs, xq, w, K, d)
        Ap = np.dot(A, p).tolist()
        r = []
        for j in range(len(Ap)):
            r.append(abs(Ap[j]-b[j]))
        y.append(max(r))
    return x, y

#Takes as input x, a list of x-values, and y, which is a list containing multiple lists yi that contain y-values corresponding
#to the x-values, and x to yi for each yi in y.
def Plot_func(x, y):
    for yi in y:
        plt.plot(x, yi)
    plt.yscale('log')
    plt.show()

#Used in problem 5 - takes as input minimum and maximum values for Nc given Nc=Ns. Iterates through these and for each value
#Solves the linear system Ap=b to find p, and calculates the inf-norm of p minus the analytical solution.
#Returns x, y where x is the Nq-values, and y contains the corresponding values of the inf-norm.
def Gen_Error_p(start, end, K, F, method, d):
    x = []
    y = []
    for i in range(start, end+1):
        print(i)
        xc = Chebyshev(i)
        b = fredholm_rhs(xc, F, d)
        xq, w = method(i**2)
        A = fredholm_lhs(xc, xc, xq, w, K, d)
        p = np.linalg.solve(A, b)
        p2 = Density(xc)
        r = []
        for j in range(len(p)):
            r.append(abs(p[j]-p2[j]))
        x.append(i)
        y.append(max(r))
    return x, y

#Generates b and the perturbed version of b at a given depth d, and plots the both in the same graph.
def Gen_Perturbed(n, F, d):
    x = []
    y = []
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F, d)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    plt.plot(xc, b, label = 'Not perturbed')
    plt.plot(xc, b2, label = 'Perturbed')
    plt.legend()
    plt.show()

#Given a depth d, and Nc = n this method solves the system Ap = b for b not perturbed and b perturbed and plots these 2
#together with the analytical solution of p in the same plot.
def Gen_plot_perturbed(n, F, K, method, d):
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F, d)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    xq, w = method(10*n)
    A = fredholm_lhs(xc, xc, xq, w, K, d)
    p1 = np.linalg.solve(A, b)
    p2 = np.linalg.solve(A, b2)
    p3 = Density(xc)
    plt.plot(xc, p1, label = 'Not perturbed')
    plt.plot(xc, p2, label = 'Perturbed')
    plt.plot(xc, p3, label = 'Analytical solution')
    plt.legend()
    plt.show()

#This method calculates a perturbed version of b, and then solves the system (ATA + Lambda*I) * p = ATb
#for Lambda = 10^i for i=-14, -13, ..., 0, 1 and for each lambda plots the numerical estimate and the analytical solution
#in the same graph. Afterwards it also plots the error as a function of lambda in a loglog plot to use for trying to find the
#optimal lambda given this specific perturbation.
def Tikhonov(n, F, K, method, d):
    xc = Chebyshev(n)
    b = fredholm_rhs(xc, F, d)
    xq, w = method(10*n)
    A = fredholm_lhs(xc, xc, xq, w, K, d)
    p1 = Density(xc)
    b2 = [x*(1+np.random.uniform(-10**-3, 10**-3)) for x in b]
    e = []
    l = []
    for i in range(-14, 2):
        print(i)
        lhs = np.dot(A.T, A) + np.dot(10**i, np.identity(n))
        rhs = np.dot(A.T, b2)
        p = np.linalg.solve(lhs, rhs)
        plt.plot(xc, p, label = 'Numerical estimate')
        plt.plot(xc, p1, label = 'Analytical solution')
        plt.legend()
        plt.show()
        r = []
        for j in range(len(p)):
            r.append(abs(p[j]-p1[j]))
        e.append(max(r))
        l.append(10**i)
    plt.plot(l, e, label = 'Error as function of lambda, d=2.5')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

#This method takes as input a file with values of xc and F(xc) and for lambda = 10^i for i = -14, -13, ..., 0, 1 and plots
#the graph given this lambda. Unfortunately we don't here have an analytical solution to plot against, so we will have to guess
#what the optimal plot looks like. Generally this was done by observing when the plots were very similar for multiple different
#lambda values in a row, as it seems as though the graph is stable for a while before and after the "optimal" lambda value usually.
#Generally this tended to be around lambda = 10^-5 or 10^-4.
def Reconstruct_Density(file, K):
    f = open(file, 'rb')
    npzfile = np.load(f)
    #sinus = lambda x: sin(5*pi*x)
    #y = [sinus(x) for x in npzfile['xc']]
    xs = Chebyshev2(len(npzfile['xc']), npzfile['a'], npzfile['b'])
    xq, w = Legendre2(len(npzfile['xc']**2), npzfile['a'], npzfile['b'])
    A = fredholm_lhs(npzfile['xc'], xs, xq, w, K, npzfile['d'])
    r = [npzfile['xc']]
    for i in range(-14, 2):
        print(i)
        lhs = np.dot(A.T, A) + np.dot(10**i, np.identity(len(npzfile['xc'])))
        rhs = np.dot(A.T, npzfile['F'])
        p = np.linalg.solve(lhs, rhs)
        plt.plot(npzfile['xc'], p)
        r.append(p)
        #plt.plot(npzfile['xc'], y, label = 'The function y=sin(5*Pi*x)')
        plt.legend()
        plt.show()
    return r

#This code snippet runs the code necessary for problem 3 and 4 (Line with Legendre in it is problem 4 and line with Trapezoid is problem 3).
'''a = Chebyshev(40)
p = Density(a)
K = lambda x, y, d: d * (d**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
x, y1 = Gen_Error(p, a, a, K, F, Legendre, 0.025)
x, y2 = Gen_Error(p, a, a, K, F, Trapezoid, 0.025)
y = [y1, y2]
Plot_func(x, y)'''

#This code snippet runs problem 5 and plots the 3 graphs at different depths as a function of Nc.
'''K = lambda x, y, d: d * (d**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
x, y1 = Gen_Error_p(5, 30, K, F, Legendre, 0.025)
x, y2 = Gen_Error_p(5, 30, K, F, Legendre, 0.25)
x, y3 = Gen_Error_p(5, 30, K, F, Legendre, 2.5)
y = [y1, y2, y3]
Plot_func(x, y)'''

#This code snippet generates an exact and perturbed graph for problem 6 (that tended to just be basically the same graph from
#the perspective of a viewer as the perturbations of 0.1% or less weren't enough to make for very visible differences.
'''F = pickle.load( open( "F.pkl", "rb" ) )
Gen_Perturbed(30, F, 0.025)'''

#This code snippet is used for problem 6 and 7 - uncomment gen_plot_perturbed for problem 6 and tikhonov for problem 7.
#Input desired d instead of 2.5 if you want some other depth.
'''K = lambda x, y, d:  d * (d**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
#Gen_plot_perturbed(30, F, K, Legendre, 2.5)
#Tikhonov(30, F, K, Legendre, 2.5)'''

#The first 4 lines will plot a bunch of graphs and such for each of the 3 files separately to try to find the optimal lambdas.
#The last half will iterate through the different plots with the values from each of the 3 files to have all the graphs in the
#same plot. In each case it will start at lambda = 10^-14 and end at lambda =10^1, printing the exponent i at each step.
'''K = lambda x, y, d:  d * (d**2 + (y-x)**2)**(-3/2)
r1 = Reconstruct_Density('q8_1.npz', K)
r2 = Reconstruct_Density('q8_2.npz', K)
r3 = Reconstruct_Density('q8_3.npz', K)

for i in range(1, len(r1)):
    print(-15+i)
    plt.plot(r1[0], r1[i], label = 'First file')
    plt.plot(r2[0], r2[i], label = 'Second file')
    plt.plot(r3[0], r3[i], label = 'Third file')
    plt.legend()
    plt.show()'''