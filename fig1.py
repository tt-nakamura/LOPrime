import numpy as np
import matplotlib.pyplot as plt
from mpmath import zeta

Z = np.vectorize(lambda z: complex(zeta(z)))
s = 1.5
x = 100
T = 50
n = 256

t = np.linspace(0, T, n+1)
z = s + t*1j
F = np.exp(z**2/2/x)
w = x**z/z * np.log(Z(z))
w1 = np.real(w)
w2 = np.real(w*F)
plt.semilogy(t,np.abs(w1), label=r'$F=1$')
plt.semilogy(t,np.abs(w2), label=r'$F=e^{z^2/(2x)}$')
plt.xlabel(r'$t = {\rm Im}\,z$', fontsize=14)
plt.ylabel(r'| Re $(x^z/z)\,\log\,\zeta(z)\,F$ |', fontsize=14)
plt.text(0, 2e-5, r'Re $z = 3/2$', fontsize=14)
plt.text(0, 1e-4, r'$x = %d$'%x, fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()


