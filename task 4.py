import numpy as np
from math import *
from matplotlib import pyplot as plt
import pickle
import scipy.special as ss
import sympy as sy

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
                A[i][j] += K(xc[i], (xq[k]+1)/2) * .5*w[k] * Lagrange_Basis(j, (xq[k]+1)/2, xs, Ns)
    return A

def Chebyshev(n):
    a = []
    for i in range(1, n+1):
        a.append(0.5 + 0.5*cos(pi*(2*i-1)/(2*n)))
    return a

"""
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
"""

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

def Gen_Error(n, p, xc, xs, K, F):
    x = []
    y = []
    b = fredholm_rhs(xc, F)
    for i in range(1, 30):
        xq, w = np.polynomial.legendre.leggauss(i)
        x.append(i)
        A = fredholm_lhs(xc, xs, xq, w, K)
        Ap = np.dot(A, p).tolist()
        r = []
        for j in range(len(Ap)):
            r.append(abs(Ap[j]-b[j]))
        y.append(max(r))
    return x, y


def Plot_func(x, y):
    plt.plot(x, y)
    plt.show()

class analytical_solution:
    """
    Evaluate the approximation to the force measurement
    corresponding to rho(x) = sin(omega*x) exp(gamma x)
    based on Nmax Taylor series terms.
    We do Taylor series expansion at (a+b)/2.
    Distance from the measurements to the mass density is d.
    """
    def __init__(self,a,b,omega,gamma,Nmax):
        import time
        """
        Initialize the object: do the Taylor series expansion etc
        """
        self.a = a
        self.b = b

        # define symbols
        x,y,u = sy.symbols('x y u', real=True)
        # define a density function we want to integrate
        rho = sy.sin(omega*y) * sy.exp(gamma*y)
        #
        # make a Taylor series expansion of this density
        # up to Nmax terms
        rho_taylor = sy.series(rho,y,(a+b)/2,Nmax).removeO()
        # Now we substitute y=u+x and represent the result as a polynomial wrt u
        pu_coeffs = rho_taylor.subs(y,u+x).as_poly(u).all_coeffs()
        self.pu_coeffs_str = []
        for coeff in pu_coeffs:
            self.pu_coeffs_str.append(str(coeff))
        # for evaluation we would like to convert these functions
        # to lambda-function, but those cannot be stored (pickled)
        # we will store the lambda functions in the following list:
        self.cns = []

    def perform_lambdification(self):
        """
        Convert the extracted Taylor series coefficients
        to efficiently evaluatable functions
        """
        x = sy.symbols('x')
        self.cns = []
        for n in range(len(self.pu_coeffs_str)):
            # extract the polynomial coefficient corresponding to u^n as a function of x
            pu_coeff_n = sy.sympify(self.pu_coeffs_str[-1-n])
            cn = sy.lambdify(x,pu_coeff_n,"numpy")
            self.cns.append(cn)

    def antideriv(selv,u,d,n):
        """
        Antiderivative of  u^n/(d^2+u^2)^1.5
        """
        return u**(n+1)  * ss.gamma(n/2+0.5)/      \
               (2 * d**3 * ss.gamma(n/2+1.5)) *    \
               ss.hyp2f1(1.5,n/2+0.5,n/2+1.5,-(u/d)**2)

    def __call__(self,x_eval,d):
        """
        Evaluate the initialized object at x_eval
        """
        if self.cns == []:
            self.perform_lambdification()
        if np.isscalar(x_eval):
            x_eval = np.array([x_eval])
        F_eval = np.zeros_like(x_eval)
        for n in range(len(self.cns)):
            F_eval = F_eval + d*self.cns[n](x_eval) * \
                (self.antideriv(self.b-x_eval,d,n)-self.antideriv(self.a-x_eval,d,n))
        return F_eval



a = Chebyshev(40)
p = Density(a)
K = lambda x, y: 0.025 * (0.025**2 + (y-x)**2)**(-3/2)
F = pickle.load( open( "F.pkl", "rb" ) )
x, y = Gen_Error(100, p, a, a, K, F)
Plot_func(x, y)