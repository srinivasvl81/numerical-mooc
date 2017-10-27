# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 07:28:38 2017
# Traffic flow assignment

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
L = 11          #km
Vmax = 80       # km/hr
rho_max = 250   # cars/km
nx = 51
nt = 20
dt = 0.001      #hrs

x = np.linspace(0, L, nx)

#c =1        # Wave speed

# Initial conditions
rho0 = np.ones(nx)*10
rho0[10:20] = 50

#rho = np.asanyarray(rho0)


# Plotting the wave
pp.figure(figsize=(10, 6))
pp.title('Traffic flow assignment', fontsize=18)
pp.grid(True)
pp.plot(x, rho0, color = '#003366', ls ='--', lw=3)
pp.xlim(0,11)
pp.ylim(0, 60);
#
# Iteration
for n in range(1, nt):      # Time loop
    rho_n = rho.copy()
    """ At every time step, rho_n becomes the previous time step
    values. For every time step, space loop runs and determines
    next time step values.
    """           
#    u[1:] = un[1:] - un[1:]*dt/dx*(un[1:] - un[0:-1]) #in one swoop   
    for i in range(1, nx):      # Spatial loop
        # FTBS eq
        rho[i] = rho_n[i] - Vmax*dt/dx* (1 - 2*rho_n[i]/rho_max)*(rho_n[i]-rho_n[i-1])
    
    u[0] = 1.0
    
## Plotting
#pp.plot(x, u, color = '#003366', ls ='--', lw=3)
#pp.ylim(0, 2.5);
#
#elapsed = time.time() - start 
#print ("The elapsed time is: ", elapsed*1000.,"milli-seconds")