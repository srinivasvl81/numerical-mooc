# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 10:01:17 2017

@author: srinivas
"""

from math import sin, cos, log
import numpy as np
from matplotlib import pyplot as pp
import time
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start_time = time.time()

# model parameters
g = 9.8     # gravity in m/s
v_t = 4.9   # trim velocity in m/s
C_D = 1/5.0     # Drag coefficient 
C_L = 1.0       # Lift coefficient

# Setting initial conditions
v0 = 6.5    
theta0 = -0.1       # initial angle of the trajectory
x0 = 0.0            # Horizantal position is arbitary
y0 = 2.0            # initial altitude

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

# Euler step 
def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u)

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

# set time discretization
T = 15.0            # final time
dt = 0.01           # set time-increment
N = int(T*1.0/dt) + 1   # Number of time-steps

# set the initial conditions
u_euler = np.empty((N, 4))
u_rk2 = np.empty((N, 4))

# initialize the array containing the solution for each time-step
u_euler[0] = np.array([ v0, theta0, x0, y0])
u_rk2[0] = np.array([ v0, theta0, x0, y0])

# use a loop to call rk2_step()
for n in range(N-1):
    u_euler[n+1] = euler_step(u_euler[n], f, dt)
    u_rk2[n+1] = rk2_step(u_rk2[n], f, dt)
    
# assigning the position of glider in time to separate arrays
x_euler = u_euler[:, 2]
y_euler = u_euler[:, 3]
x_rk2 = u_rk2[:,2]
y_rk2 = u_rk2[:,3]

# Get the index of the element where y- coordinate i.e. altitude becomes negative
idx_negative_euler = np.where(y_euler<0.0)[0]
if len(idx_negative_euler)==0:
    idx_ground_euler = N-1
    print ('Euler integration has not touched ground yet!')
else:
    idx_ground_euler = idx_negative_euler[0]
    
idx_negative_rk2 = np.where(y_rk2<0.0)[0]
if len(idx_negative_rk2)==0:
    idx_ground_rk2 = N-1
    print ('Runge-Kutta integration has not touched ground yet!')
else:
    idx_ground_rk2 = idx_negative_rk2[0]
    
# Checking if the paths match
print 'Are the x-values close ? {}'.format(np.allclose(x_euler, x_rk2))
print 'Are the y-values close ? {}'.format(np.allclose(y_euler, y_rk2))

# plot the glider path
pp.figure(figsize=(10,6))
pp.subplot(121)
pp.grid(True)
pp.xlabel('$x$')
pp.ylabel('$y$')
pp.plot(x_euler[:idx_ground_euler], y_euler[:idx_ground_euler], 'k-', label='Euler')
pp.plot(x_rk2[:idx_ground_rk2], y_rk2[:idx_ground_rk2], 'r--', label='RK-2')
pp.title('Distance Travelled: {:.3f}'.format(x_rk2[idx_ground_rk2-1]))
pp.legend();

# for a closer look!
pp.subplot(122)
pp.grid(True)
pp.xlabel('$x$')
pp.ylabel("$y$")
pp.plot(x_euler, y_euler, 'k-', label='Euler')
pp.plot(x_rk2, y_rk2, 'r--', label='Runge-Kutta-2')
pp.xlim(0,5)
pp.ylim(1.8, 2.5);

# Grid convergence
# Using a for loop to compute the solution on different grids
dt_values= np.array([0.1, 0.05, 0.01, 0.005, 0.001])

u_values = np.empty_like(dt_values, dtype=np.ndarray)

for i,dt in enumerate(dt_values):
    N = int(T*1.0/dt) + 1        # Number of time steps
    # discretizing the time t
    t = np.linspace(0.0, T, N)
    
    # Initialize the array containing the solution for each time step
    u = np.empty((N, 4))
    u[0] = np.array([v0, theta0, x0, y0])
    
    #time loop
    for n in range(N-1):
        u[n+1] = rk2_step(u[n], f, dt)
    
    # store the value of u related to one grid
    u_values[i] = u

# Compute diffgrid
diffgrid = np.empty_like(dt_values)
for i, dt in enumerate(dt_values):
    diffgrid[i] = get_diffgrid(u_values[i], u_values[-1], dt)
    

# Plot the gird convergence using matplotlib loglog()
pp.figure(figsize=(6,6))
pp.grid(True)
pp.xlabel(r'$\Delta t$', fontsize = 18)
pp.ylabel(r'$L_1$ - norm of grid differences', fontsize =18)
pp.xlim(1e-4, 1)
pp.ylim(1e-4, 1)
pp.axis('equal')
pp.loglog(dt_values[:-1], diffgrid[:-1], color ='k',ls='--',lw=2,marker='o');


# Checking the convergence rate
r = 2
h = 0.001
dt_values = np.array([h, r*h, r**2*h])
print dt_values
u_values = np.empty_like(dt_values, dtype = np.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T*1.0/dt) + 1       # Number of time-steps
    
    # Discretize the time
    t = np.linspace(0.0, T, N)
    
    # Initialize the array containing the solution for each time step
    u = np.empty((N, 4))
    u[0] = np.array([v0, theta0, x0, y0])
    
    # Time loop
    for n in range(N-1):
        # Call rk2_step() 
        u[n+1] = rk2_step(u[n], f, dt)
    
    # store the values of u related to one grid
    u_values[i] = u

# Calculating the order of convergence
alpha = (log(get_diffgrid(u_values[2], u_values[1], dt_values[2]))-\
        log(get_diffgrid(u_values[1], u_values[0], dt_values[1])))/log(r)*1.0

print 'The order of convergence is alpha = {:.3f} '.format(alpha)


elapsed_time = time.time() - start_time

print 'The time taken by Runge-Kutta 2 method is', elapsed_time,'seconds' 
    
