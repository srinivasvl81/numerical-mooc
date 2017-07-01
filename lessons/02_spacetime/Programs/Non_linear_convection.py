# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 17:21:39 2017
1D Non-linear convection
@author: srinivas
"""

# Importing modules
import numpy as np
import time
from matplotlib import pyplot as pp
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start = time.time()
# Parameters
nx = 41
dx = 2.0/(nx-1)
nt = 10
dt = .02
c =1        # Wave speed
x = np.linspace(0, 2, nx)

# Defining the bounds of the wave
u = np.ones(nx)
lbound = np.where(x >= 0.5)
ubound = np.where(x <= 1.0)

u[np.intersect1d(lbound, ubound)] = 2.0

# Plotting the wave
pp.figure(figsize=(10, 6))
pp.title('1D Non-Linear Convection', fontsize=18)
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
    u[1:] = un[1:] - un[1:]*dt/dx*(un[1:] - un[0:-1]) #in one swoop   
    #for i in range(1, nx):      # Spatial loop
        #u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1])  # FTBS eq
    u[0] = 1.0
# Plotting
pp.plot(x, u, color = '#003366', ls ='--', lw=3)
pp.ylim(0, 2.5);

elapsed = time.time() - start 
print ("The elapsed time is: ", elapsed*1000.,"milli-seconds")