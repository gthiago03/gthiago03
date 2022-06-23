import numpy as np
import scipy
import sympy

import scipy.sparse as sps

import math

import porepy as pp

###
np.random.seed(42)
base = 4
domain = np.array([1, 1])
basedim = np.array([base, base])
num_refs = 7
pert = 0.5
ref_rate=2
### End of parameter definitions

# Analytical solution
x, y = sympy.symbols('x y')

kd = 1 #diagonal permeability
kc = 0.1 #cross permeability

g1 = 10 
g2 = 1

l = 1

u0 = 1 #reference pressure

ht = domain[1]
hi = domain[1]/2

# discontinuous + smooth
a1 = 1
a2 = 0
u = u0 + a1 * sympy.Piecewise(((ht-y)*g1, y>=hi), ((ht-hi)*g1+(hi-y)*g2, y<hi)) + a2 * sympy.sin(x) * sympy.cos(y)

gx = sympy.diff(u, x)
gy = sympy.diff(u, y)

u_f = sympy.lambdify((x, y), u, 'numpy')
gx_f = sympy.lambdify((x, y), gx, 'numpy')
gy_f = sympy.lambdify((x, y), gy, 'numpy')

