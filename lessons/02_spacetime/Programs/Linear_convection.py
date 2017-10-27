# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:21:48 2017
1D linear convection
@author: srinivas
"""

# Importing modules
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# Parameters
nx = 50
dx = 0.2#2.0/(nx-1)
print dx
nt = 10
dt = .05
c =4       # Wave speed
x = np.linspace(0, 2, nx)

# Defining the bounds of the wave
u = np.ones(nx)
lbound = np.where(x >= 0.5)
ubound = np.where(x <= 1.0)

#print (lbound)
#print (ubound)


bounds = np.intersect1d(lbound, ubound)

u[bounds] = 2.0
#print (u)

# Plotting the wave
pp.figure(figsize=(10, 6))
pp.title('1D Convection', fontsize=18)
pp.grid(True)
pp.plot(x, u, color = '#003366', ls ='--', lw=3)
pp.ylim(0, 2.5);

# Iteration
for n in range(1, nt):      # Time loop
    un = u.copy()
    """ At every time step, un becomes the previous time step
    values. For every time step, space loop runs and determines
    next time step values.
    """           
    for i in range(1, nx):      # Spatial loop
        u[i] = un[i] - c*dt/dx*(un[i]-un[i-1])  # FTBS eq

# Plotting
pp.plot(x, u, color = '#003366', ls ='--', lw=3)
pp.ylim(0, 2.5);


