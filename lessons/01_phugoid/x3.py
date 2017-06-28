# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 21:25:54 2017

@author: srinivas
"""


# Importing and setting
from math import sin, cos, log, ceil
import numpy
import time
from matplotlib import pyplot
from matplotlib import rcParams 
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start_time= time.time()
# Model parameters
g = 9.8                     # gravity in m/s^2
v_t = 30.0                  # trim velocity in m/s
C_D = 1/40.0                # Drag coefficient
C_L = 1.0                   # Lift coefficient

# ---------- Initial conditions ----------
v0 = v_t                    # starting at the trim velocity
theta0 = 0.                  # Initial angle of trajectory
x0 = 0.                      # arbitary horizantal position
y0 = 1000.0                    # initial altitude

# Defining the function f() to return the phugoid RHS
def f(u):
    """ PARAMETERS:
    u = array of float
        array containing the solution at time n
        
    RETURNS
    dudt = array of float
            array containing RHS given u.
    """
    v = u[0]
    theta = u[1]
    x = u[2]
    y = u[3]
    return numpy.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                      -g*cos(theta)/v + g/v_t**2*v,
                      v*cos(theta),
                      v*sin(theta)])

# Defining a function Euler step
def euler_step(u, f, dt):
    u_star = u + 0.5*dt*f(u)
    return u + dt* f(u_star)     
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

    
# Discretizing the time
T = 100                     # final time
dt = 0.1                   # time increment
N = int(T/dt) + 1           # Number of time steps

# initializing the solution array 
u = numpy.empty((N, 4))
u[0] = numpy.array([v0, theta0, x0, y0])  # filling the first elements

# Time loop - Euler method
for n in range(N - 1):
    u[n+1] = euler_step(u[n], f, dt)

# get the glider's position with respect to the time
x = u[:, 2]
y = u[:, 3]

# visualization of the path
pyplot.figure(figsize=(10, 6))
pyplot.grid(True)
pyplot.xlabel(r'x', fontsize=18)    
pyplot.ylabel(r'y', fontsize=18)
pyplot.title('Glider Trajectory, Flight time = %.2f' % T, fontsize=18)
pyplot.plot(x, y, 'k-', lw = 2);

# Grid Convergence
dt_values = numpy.array([0.1, 0.05, 0.01, 0.005, 0.001])

u_values = numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt) + 1   # Number of time steps
    
    # initialize the array containing the solution
    u = numpy.empty((N, 4))
    u[0] = numpy.array([v0, theta0, x0, y0])
    
    # time loop
    for n in range(N-1):
        u[n + 1] = euler_step(u[n], f, dt)        #calling euler_step
        
    # store the value of u related to one grid
    u_values[i] = u

# Grid difference between coarse and fine grids
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
   
    grid_size_ratio = ceil(N_fine*1.0/N_current)
    
    diffgrid = dt * numpy.sum( numpy.abs(\
            u_current[:,2]- u_fine[::grid_size_ratio,2])) 
    
    return diffgrid

# Compute the difference between one grid solution and the finest grid
diffgrid = numpy.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    print('dt = {}'.format(dt))

    ### call the function get_diffgrid() ###
    diffgrid[i] = get_diffgrid(u_values[i], u_values[-1], dt)
    

# log-log plot of the grid differences
pyplot.figure(figsize=(6,6))
pyplot.grid(True)
pyplot.xlabel('$\Delta t$', fontsize=18)
pyplot.ylabel('$L_1$-norm of the grid differences', fontsize=18)
pyplot.axis('equal')
pyplot.loglog(dt_values[:-1], diffgrid[:-1], color='k', ls='-', lw=2, marker='o');

# Order of convergence
r = 2
h = 0.001

dt_values2 = numpy.array([h, r*h, r**2*h])

u_values2 = numpy.empty_like(dt_values2, dtype=numpy.ndarray)

diffgrid2 = numpy.empty(2)

for i, dt in enumerate(dt_values2):
    
    N = int(T*1.0/dt) + 1   # number of time-steps

    # initialize the array containing the solution for each time-step
    u = numpy.empty((N, 4))
    u[0] = numpy.array([v0, theta0, x0, y0])

    # time loop
    for n in range(N-1):

        u[n+1] = euler_step(u[n], f, dt)         ### call euler_step() ###
    
    # store the value of u related to one grid
    u_values2[i] = u
    

#calculate f2 - f1
diffgrid2[0] = get_diffgrid(u_values2[1], u_values2[0], dt_values2[1])

#calculate f3 - f2
diffgrid2[1] = get_diffgrid(u_values2[2], u_values2[1], dt_values2[2])

# calculate the order of convergence
p = (log(diffgrid2[1]) - log(diffgrid2[0])) / log(r)

print('The order of convergence is p = {:.3f}'.format(p));
elapsed_time = time.time() - start_time
print 'The elapsed time in Euler method is:', elapsed_time, 'seconds'
