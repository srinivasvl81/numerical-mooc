# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 18:54:14 2017

@author: srinivas
"""

import numpy
from matplotlib import pyplot
T=100.0
dt=0.02
N=int(T/dt)+1   # Number of grid points
t=numpy.linspace(0.0, T, N)

# initial conditions
z0=100.  # altitude
b0=10.   # upward velocity resulting from gust
zt=100.
g=9.81   #acceleration due to gravity

u=numpy.array([z0,b0])

#initialize an array to hold the changing elevation values
z=numpy.zeros(N)   
z[0]=z0

# Time looping using Euler's method
for i in range(1,N):
    u = u + dt*numpy.array([u[1], g*(1.-u[0]/zt)])
    z[i] = u[0]
 
"""
# Plotting the numerical solution
pyplot.figure(figsize=(10,4))   # set plot size
pyplot.ylim(40,160)             # y-axis plot limits
pyplot.tick_params(axis='both', labelsize=14)   # increase font size for ticks
pyplot.xlabel('t', fontsize=14) # x-label
pyplot.ylabel('z', fontsize=14) # y-label
pyplot.plot(t, z, 'k-');
"""

# Plotting the exact solution
z_exact = b0*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
            (z0-zt)*numpy.cos((g/zt)**.5*t)+zt

pyplot.figure(figsize=(10,5))
pyplot.ylim(40,160)             #y-axis plot limits
pyplot.tick_params(axis='both', labelsize=14) #increase font size for ticks
pyplot.xlabel('t', fontsize=14) #x label
pyplot.ylabel('z', fontsize=14) #y label
pyplot.plot(t,z)
pyplot.plot(t, z_exact)
pyplot.legend(['Numerical Solution','Analytical Solution']);

# -------------CONVERGENCE-------------
# Time increment array
dt_values=numpy.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0001])

# array will contain solution of each grid
z_values=numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt) + 1   # Number of time steps
    t = numpy.linspace(0.0, T, N)    # time discretization using linspace
    
    # initial conditions
    u=numpy.array([z0, b0])
    z=numpy.empty_like(t)
    z[0]=z0
    
    #time loop - Euler method
    for n in range(1, N):
        # computer next solution using Euler method
        u = u + dt*numpy.array([u[1], g*(1.-u[0]/zt)])
        z[n] = u[0] # store the elevation at time step n+1
    
    z_values[i] = z.copy()
    
# Calculating the error
def get_error(z, dt):
    """Returns the error relative to analytical solution using L-1 norm.
    
    Parameters
    ----------
    z : array of float
        numerical solution.
    dt : float
        time increment.
        
    Returns
    -------
    err : float
        L_{1} norm of the error with respect to the exact solution.
    """
    N = len(z)
    t = numpy.linspace(0.0, T, N)
    z_exact = b0*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
            (z0-zt)*numpy.cos((g/zt)**.5*t)+zt
    
    return dt*numpy.sum(numpy.abs(z-z_exact))    

# Error values
error_values=numpy.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    # Calling the function get_error#
    error_values[i] = get_error(z_values[i], dt)
    
# Plotting dt vs error
pyplot.figure(figsize=(10,8))
pyplot.tick_params(axis='both', labelsize = 14)    #increase tick font size
pyplot.grid(True)                               # Turn on the grid
pyplot.xlabel('$\Delta t$', fontsize = 16)      # x label
pyplot.ylabel('Error', fontsize = 16)           # y label
pyplot.loglog(dt_values, error_values, 'ko-')   # log-log plot
pyplot.axis('equal')                        # makes axes scale equally


# Done !


            