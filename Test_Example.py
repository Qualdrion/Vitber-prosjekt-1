import sympy as sy
import numpy as np
import scipy.special as ss
import math

"""
In this example we construct an approximation
to the integral

int_a^b K(x,y) rho(y) dy

as described in the appendix A of the project desctiption;
that is, we perform a Taylor series expansion of rho(),
and then use the special function representation of the
integral over polynomial.
"""

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


if __name__=='__main__':
    import matplotlib.pyplot as plt
    import time
    import pickle

    #
    # integration limits
    a = 0
    b = 1
    # d is the distance from the measurements in the kernel
    d = 2.5E-02
    # we use rho(y) = exp(gamma*y)*sin(omega*y) as an example
    gamma = -2
    omega = 3*np.pi
    # maximal order of the Taylor series expansion
    Nmax  = 75
    # evaluate the integral expression at x_eval
    N_eval = 50
    x_eval = np.linspace(a,b,N_eval)
    start_t= time.time()
    F      = analytical_solution(a,b,omega,gamma,Nmax)
    end_t  = time.time()
    print("Initialization took %f s." % (end_t-start_t))
    # Since the initialization is slow, we can store
    # the object in a file (pickle it).
    pickle.dump( F, open( "F.pkl", "wb" ) )

    # ... much later, in a galaxy far, far away

    start_t= time.time()
    # Now, whenever we want, we can simply load the object from file
    # and use it
    try:
        F = pickle.load( open( "F.pkl", "rb" ) )
    except:
        F = analytical_solution(a,b,omega,gamma,Nmax)

    F_eval = F(x_eval,d)
    end_t  = time.time()
    print("First evaluation took %f s" % (end_t-start_t))

    start_t= time.time()
    F_eval2= F(x_eval,d)
    end_t  = time.time()
    print("Second evaluation took %f s" % (end_t-start_t))

    # finally, plot the function
    plt.plot(x_eval,F_eval2, 'k')
    plt.xlabel('x')
    plt.ylabel('F(x)')
    plt.grid(True)
    plt.show()