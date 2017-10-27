# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 17:51:18 2017
1D Diffusion
@author: srinivas
"""
# Importing modules
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import rcParams
import time
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start = time.time()

nx = 41         # Number of spatial grid points
dx = 2.0/(nx-1)
nt = 50
nu = 0.3        # Value of viscosity
sigma = 0.2     # CFL number
dt = sigma*dx**2/nu     # time step size

x = np.linspace(0.0, 2, nx)  # spatial grid array
lbound = np.where(x >= 0.5)
ubound = np.where(x <= 1)

# Defining the HAT function
u = np.ones(nx)
u[np.intersect1d(lbound, ubound)] = 2

un = np.ones(nx)

# Loop to calculate the solution
for n in range(nt):
    un = u.copy()
    u[1:-1] = un[1:-1] +(nu*dt/dx**2)*(un[2:] + un[0:-2] - 2*un[1:-1])

# Plotting   
pp.grid(True) 
pp.plot(x, u, color='#003366', ls='--', lw=3)
pp.ylim(0, 2.5);

# Animating the plot
from matplotlib import animation as ani

fig = pp.figure(figsize = (8, 5))
ax = pp.axes(xlim = (0,2), ylim=(1,2.5))
line = ax.plot([], [], color='#003366', ls='--', lw=3)[0]
ax.set_title('1D Diffusion Equation-Solution')

# Diffusion function
def diffusion(i):
    line.set_data(x, u)
    
    un = u.copy()
    u[1:-1] = un[1:-1] + (nu*dt/dx**2)*(un[2:] + un[0:-2] - 2*un[1:-1])

# Animation 
anim = ani.FuncAnimation(fig, diffusion, frames=nt, interval=100)

#anim.save('image.mp4', fps=20, writer="avconv", codec="libx264")
anim.save('Diffusion.mp4', writer='avconv')

elapsed = time.time() - start
print('Time taken:', elapsed,'sec') 