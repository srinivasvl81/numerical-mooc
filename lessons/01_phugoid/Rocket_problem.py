# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:47:32 2017

@author: srinivas
"""

# Coding assignment: Rocket flight

# Importing modules
from math import log, pi
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# Rocket data
ms = 50.0         # mass of rocket shell in kgs
g = 9.81        # acceleration due to gravity in m/s^2
rho = 1.091     # density in kg/m^3
r = 0.5         # rocket cross section radius in m
A = pi*r**2     # Rocket cross section area in m^2/s
ve = 325.0        # Exit velocity of mass in m/s
C_D = 0.15      # Drag coefficient
mp0 = 100.0       # Initial mass of propellent in kgs at time T = 0
D = 0.5*rho*A*C_D

# initial conditions
t0 = 0.      # initial time
h0 = 0.       # initial height = ground level
v0 = 0.       # initial velocity; at rest

# Defining the function f() = RHS of rocket equations
def f(u):
    """ Parameters:
    u = array of float, array containing the solution at time n
    
    Returns:
    du/dt = array of float; array containing RHS of given u
    """
    t = u[0]
    h = u[1]
    v = u[2]
    if t<=5.0:
        return np.array([1.0, v*1.0, (-g-(D*v*abs(v)-20.0*ve)/(150.0-20.0*t))])
    elif t>5.0:
        return np.array([1.0, v*1.0, (-g - (D/ms)*v*abs(v) )])
    
# Defining a function Euler step
def euler_step(u, f, dt):
    return u + dt * f(u)     
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equations.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """

# Discretizing the time
T = 40                   # final time in s
dt = 0.1                   # time increment in s
N = int(T*1.0/dt) + 1           # Number of time steps

# initializing the solution array 
u = np.empty((N, 3))
u[0] = np.array([t0, h0, v0])  # filling the first elements

#Storing the results to a text file
res = open('results.txt', 'w')
line = 'Time-----Height-----Velocity'  
res.write(line)
res.write("\n%.3f" % u[0,0])
res.write("\t%.3f" % u[0,1])
res.write("\t%.3f" % u[0,2])
#print line
#print u[0]

# Time loop - Euler method
for n in range(N-1):
    # Writing to a text file at each and every iteration
    u[n+1] = euler_step(u[n], f, dt)
    res.write("\n%.3f" % u[n+1,0])
    res.write("\t%.3f" % u[n+1,1])
    res.write("\t%.3f" % u[n+1,2])
    #print u[n+1]

    
res.close()
# get the rocket's position with respect to the time
t = u[:, 0]
h = u[:, 1]


# visualization of the path - plotting
pp.figure(figsize=(10, 6))
pp.grid(True)
pp.xlabel(r't', fontsize=18)    
pp.ylabel(r'h', fontsize=18)
pp.title('Rocket Trajectory, Flight time = %.2f' % T, fontsize=18)
pp.plot(t, h, 'k-', lw = 2);   
    
        