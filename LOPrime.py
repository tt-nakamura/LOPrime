# references:
#   J. C. Lagarians and A. M Odlyzko
#    "Computing pi(x): Analytic Method"
#     Journal of Algorithms 8 (1987) 173
#   R. Crandall and C. Pomerance
#    "Prime Numbers: A Computational Perspective"
#     section 3.7.2

import numpy as np
from mpmath import zeta
from sympy import sieve
from scipy.integrate import quad
from scipy.special import erfc

s = 1.5 # real part of integral path

def LOPrime_(x,T,**kw):
    """ count number of primes in [1,x]
        by Lagarias-Odlyzko method
    T = upper limit of integral
    kw = keyword arguments passed to scipy.integrate.quad
    """
    def f(t):
        z = s + 1j*t
        Z = complex(zeta(z)) # Riemann zeta function
        return np.real(x**z/z * np.log(Z))

    I = quad(f,0,T,**kw)[0]/np.pi
    m = np.arange(2, x.bit_length())
    c = [len(list(sieve.primerange(1,p))) for p in x**(1/m)]
    return I - np.sum(c/m)

K = 3 # width of smoothing function

def LOPrime(x,T,**kw):
    """ count number of primes in [1,x]
        by Lagarias-Odlyzko method with gaussian
        smoothing function for convergence factor
    T = upper limit of integral
    kw = keyword arguments passed to scipy.integrate.quad
    """
    def f(t):
        z = s + 1j*t
        Z = complex(zeta(z)) # Riemann zeta function
        F = np.exp(z**2/2/x) # convergence factor
        return np.real(x**z/z * F * np.log(Z))

    I = quad(f,0,T,**kw)[0]/np.pi
    a = np.sqrt(2/x)
    for m in range(1, x.bit_length() + 1):
        p0 = x**(1/m)
        k = np.exp(K*a/m)
        p = np.fromiter(sieve.primerange(p0/k, p0*k), np.int)
        v = m*np.log(p/p0)/a
        I += np.sum((p<p0) - erfc(v)/2)/m
        if m>1:
          I -= len(list(sieve.primerange(1,p0)))/m
    return I
