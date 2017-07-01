# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 21:19:33 2017
CFL condition
@author: srinivas
"""

# Importing modules
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# Defining a function for linear convection
def linear_conv(nx):
    """ Solve the 1D linear convection equation
    * wave speed c = 1
    * domain limits are x E [0, 2]
    * 20 time steps are assumed
    * dt = 0.025
    
    Produces a plot of the results
    
    Parameters: nx = no. of internal grid points
    
    Returns: None
    """
    dx = 2.0/(nx-1)
    nt = 20
    c = 1.0
    sigma = .8  # CFL number defined
    dt = sigma*dx/c      # time step determined from CFL no.
    x= np.linspace(0.0, 2, nx)
    u = np.ones(nx)
    
    lbound = np.where(x >= 0.5)
    ubound = np.where(x <= 1)
    u[np.intersect1d(lbound, ubound)] = 2
    
    un = np.ones(nx)    
    
    for n in range(nt):
        un = u.copy()
        u[1:] =  un[1:] - c*dt/dx*(un[1:] - un[0:-1])
        u[0] = 1.0
    
    pp.grid(True)    
    pp.plot(x, u, color = 'r', ls='--', lw = 2)
    pp.ylim(0, 2.5)

linear_conv(500)    
        
