import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from sympy import primerange

x = 1000
marker = ['.','x','+']
k = 0.2
x1,x2 = x*(1-k),x*(1+k)
a = np.sqrt(2/x)
i = 0
for m in range(1, x.bit_length()+1):
    p = primerange(x1**(1/m), x2**(1/m))
    u = np.asarray([x**m for x in p])
    if len(u)==0: continue;
    c = erfc(np.log(u/x)/a)/2
    y = ((u<x) - c)/m
    plt.plot(u,y,marker[i],label='m=%d'%m)
    i += 1

plt.ylim([-0.5, 0.5])
plt.xlabel(r'$p^m$ = (prime number)$^m$', fontsize=14)
plt.ylabel(r'$\{\theta(x-p^m) - c(p^m,x)\}/m$', fontsize=14)
plt.legend(fontsize=14)
plt.text(x1, -0.45, r'$x = %d$'%x, fontsize=14)
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
