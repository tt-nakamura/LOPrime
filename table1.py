import numpy as np
#from LOPrime import LOPrime_ as pi # table1
from LOPrime import LOPrime as pi # table2

for T in range(10,121,10):
    print(T,
          pi(100, T, epsabs=1e-2),
          pi(1000,T, epsabs=1e-2))
