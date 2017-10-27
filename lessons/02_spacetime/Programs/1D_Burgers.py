
# author @ VL Srinivas
# 1D Burgers equation

# Importing modules
import numpy as np
import sympy as sp
import time
from matplotlib import pyplot as pp
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

start= time.time()
#from sympy import init_printing
#init_printing()

x, nu, t = sp.symbols('x, nu, t')
phi = sp.exp(-(x-4.0*t)**2/(4.0*nu*(t+1))) + sp.exp(-(x-4.0*t-2.0*np.pi)**2.0/(4.0*nu*(t+1)))

phi_prime = phi.diff(x)
# print phi_prime

from sympy.utilities.lambdify import lambdify

#initial condition
u = -2.0*nu*(phi_prime*1.0/phi) + 4.0
# print u

u_lamb = lambdify((t, x, nu), u)        # puts the values of t, x, nu in u and solves 
#print("The value of u at t=1, x=4, nu=3 is {}.".format(u_lamb(1,4,3)))

# Variable declaration
nx = 101       # no. of spatial grid points
nt = 100        # no. of temporal grid points
dx = 2.0*np.pi/(nx-1)     # spatial step size
nu = 0.07        # diffusion coefficient
sigma = 0.1     # CFL number
dt = sigma*dx**2/nu     # time step size

x = np.linspace(0.0, 2.0*np.pi, nx)         # array of spatial grid points
un = np.empty(nx)
t = 0

u = np.asarray([u_lamb(t, x0, nu) for x0 in x])


# Plotting
pp.figure(figsize=(8,5), dpi=100)
pp.plot(x, u, color = '#003366', ls ='--', lw=3)
pp.xlim([0, 2.0*np.pi])
pp.ylim([0, 10])
pp.xlabel('x')
pp.ylabel('velocity-u')
pp.title('1D Burgers Equation - Initial Condition')

# Applying FTBS scheme with periodic BC
for n in range(nt):     # time loop
    un = u.copy()
    
    for i in range(nx-1):       # space loop
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1]) + \
            nu*dt/dx**2*(un[i+1] -2*un[i] + un[i-1])
            
    # periodicity
    u[-1] = un[-1] - un[-1]*dt/dx*(un[-1]-un[-2]) + \
            nu*dt/dx**2*(un[0] -2*un[-1] + un[-2])

# Analytical solution
u_analytical = np.asarray([u_lamb(nt*dt, xi, nu) for xi in x])

# Plotting
pp.figure(figsize=(8,5), dpi=100)
pp.plot(x, u, color = '#003366', ls ='--', lw=3, label='Computational')
pp.plot(x, u_analytical, label='Analytical')
pp.xlim([0, 2.0*np.pi])
pp.ylim([0, 10])
pp.grid(True)    
pp.legend();

# Animation 
from matplotlib import animation as atm

u = np.asarray([u_lamb(t, x0, nu) for x0 in x])

fig = pp.figure(figsize=(8,6))
pp.title('1D Burgers Equation')
ax = pp.axes(xlim=(0,2.0*np.pi), ylim=(0, 10))
line = ax.plot([], [], color='#003366',ls='--', lw=2)[0]
line2 = ax.plot([], [], 'k-', lw=2)[0]
ax.legend(['Numerical', 'Analytical'])
ax.set_title('1D Burgers Equation-Solution')


# Defining the 1D burgers function 
def burgers(n):
    un = u.copy()
    
    for i in range(nx-1):       # space loop
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1]) + \
            nu*dt/dx**2*(un[i+1] -2*un[i] + un[i-1])
            
    # periodicity
    u[-1] = un[-1] - un[-1]*dt/dx*(un[-1]-un[-2]) + \
            nu*dt/dx**2*(un[0] -2*un[-1] + un[-2])
    
    u_analytical = np.asarray([u_lamb(n*dt, xi, nu) for xi in x])
    line.set_data(x,u)
    line2.set_data(x, u_analytical)

anim = atm.FuncAnimation(fig, burgers, frames=nt, interval=100)

anim.save('Burgers_1D.mp4', writer='avconv')  

elapsed = time.time() - start
print('Time taken:', elapsed,'sec')   











