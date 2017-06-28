# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:29:09 2017

@author: srinivas
"""
# MULTI-STEP METHOD USING LEAP FROG SCHEME FOR PHUGOID MOTION PROBLEM
# importing all modules
from math import sin, cos, log
import numpy as np
from matplotlib import pyplot as pp
import time
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start_time = time.time()

# Leap frog function
def leapfrog_step(unm1, u, f, dt):    
    """Returns the solution time-step n+1) using Euler's method.
        
        Parameters
        ----------
        unm1 : array of float
            solution at time-step n-1.
        u : array of float
            solution at time-step n.walked 
        f : function
            function to compute the right hand-side of the system of equation.
        dt : float
            time-increment.
        
        Returns
        -------# plot the gliders path
pp.figure(figsize=(11,8))
pp.subplot(121)
pp.grid(True)
pp.xlabel('$x$')
pp.ylabel('$y$')
pp.plot(x_leapfrog[: idx_ground_leapfrog], y_leapfrog[:idx_ground_leapfrog], color='k',ls=':', lw=2)
pp.xlim(0,5)
pp.ylim(1.8, 2.5);

        u_n_plus_1 : array of float
            solution at time-step n+1.
    """
    return unm1 + 2*dt*f(u)

# Modified Euler method or 2nd order Runge-Kutta method
def rk2_step(u, f, dt):
    """ Returns the solutions at the next time step using 2nd order Runge-Kutta method
    Parameters:
    u : array of float; solution at the previous time step
    f : function to compute the RHS of the phugoid system of equations
    dt : float; time-increment
    
    Returns:
    u(n+1) : array of float; solution at the next time step
    """
    u_star = u + 0.5*dt*f(u)
    return u + dt* f(u_star)

# RHS of phugoid system of equations
def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """    
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return np.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])

# Differences with respect to fine grid
def get_diffgrid(u_current, u_fine, dt):
    """Returns the difference between one grid and the fine one using L-1 norm.
    
    Parameters
    ----------
    u_current : array of float
        solution on the current grid.
    u_finest : array of float
        solution on the fine grid.
    dt : float
        time-increment on the current grid.
    
    Returns
    -------
    diffgrid : float
        difference computed in the L-1 norm.
    """
    
    N_current = len(u_current[:,0])
    N_fine = len(u_fine[:,0])
   
    grid_size_ratio = int(np.ceil(N_fine*1.0/N_current))
    
    diffgrid = dt * np.sum( np.abs(\
            u_current[:,2]- u_fine[::grid_size_ratio,2])) 
    
    return diffgrid
    
# model parameters:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9    # trim velocity in m s^{-1}   
C_D = 1/5.0  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1

### set initial conditions ###
v0 = 6.5     # start at the trim velocity (or add a delta)
theta0 = -0.1 # initial angle of trajectory
x0 = 0.0     # horizotal position is arbitrary
y0 = 25.0     # initial altitude

# set time-increment and discretize the time
T  = 36.0                           # final time
dt = 0.01                             # set time-increment
N  = int(T*1.0/dt) + 1                   # number of time-steps

# set the initial conditions
u_leapfrog = np.empty((N, 4))

# initialize the array containing the solution for each time-step
u_leapfrog[0] = np.array([v0, theta0, x0, y0])

# first step using Runge-Kutta 2 method
u_leapfrog[1] = rk2_step(u_leapfrog[0], f, dt)

# Using for-loop to call leapfrog_step function
for n in range(1, N-1):
    u_leapfrog[n+1] = leapfrog_step(u_leapfrog[n-1], u_leapfrog[n], f, dt)
    
# get the glider position in time
x_leapfrog = u_leapfrog[:,2]
y_leapfrog = u_leapfrog[:,3]

# get the index of element y where the altitude becomes negative
idx_negative_leapfrog = np.where(y_leapfrog < 0)[0]

if len(idx_negative_leapfrog) == 0:
    idx_ground_leapfrog = N -1
    print 'The glider has not reached the ground yet!'
else:
    idx_ground_leapfrog = idx_negative_leapfrog[0]

# plot the glider path
pp.figure(figsize=(11,8))
pp.subplot(121)
pp.grid(True)
pp.xlabel('$x$')
pp.ylabel('$y$')
pp.plot(x_leapfrog[:idx_ground_leapfrog], y_leapfrog[:idx_ground_leapfrog], color='k', ls='-', lw=2)
pp.title('distance traveled: {:.3f}'.format(x_leapfrog[idx_ground_leapfrog-1]), fontsize=18);

# Let's take a closer look!
pp.subplot(122)
pp.grid(True)
pp.xlabel('$x$')
pp.ylabel('$y$')
pp.plot(x_leapfrog[:idx_ground_leapfrog], y_leapfrog[:idx_ground_leapfrog], color='k', ls=':', lw=2)
pp.xlim(115, 125)
pp.ylim(1, 2);
  
# Check convergence rate
r = 2
h = 0.001

dt_values = np.array([h, r*h, r**2*h])

u_values = np.empty_like(dt_values, dtype=np.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T*1.0/dt) + 1   # no of time steps
    
    ## discretize the tine
    t = np.linspace(0.0, T, N)
    
    # initialize the array containing the solution for each time step
    u = np.empty((N, 4))
    u[0] = np.array([v0, theta0, x0, y0])
    
    #time loop 
    u[1] = rk2_step(u[0], f, dt)
    for n in range(1, N-1):
        u[n+1] = leapfrog_step(u[n-1], u[n], f, dt)
    
    # store the values of u related to one grid
    u_values[i] = u

# calculate the order of convergence
alpha = (log(get_diffgrid(u_values[2], u_values[1], dt_values[2]))\
        -log(get_diffgrid(u_values[1], u_values[0], dt_values[1])))/log(r)

print 'The order of convergence is alpha ={:.3f}'.format(alpha)

elapsed_time = time.time() - start_time

print 'The time taken by Multistep method is', elapsed_time,'seconds'      
    
