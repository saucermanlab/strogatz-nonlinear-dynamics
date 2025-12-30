# strogatz_ch6: phase plane
import numpy as np
import matplotlib.pyplot as plt
#import scipy.integrate
#import scipy.optimize
#import sympy as sp

#%% Fig 6.1.4 saddle node system
def xdot(t, x):
    return np.array([x[0]+np.exp(-x[1]),-x[1]])

# create grid for phase portrait
x_vals = np.linspace(-2, 2, 100)
y_vals = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

# nullclines, solved symbolically in wolframalpha
# ** nullclineX seems to be a bit off **
nullclineX = -np.log(-x_vals) # 0 - x+exp(-y), solved for y
nullclineY = x_vals*0          # 0 + y, solved for y

# plotting
plt.streamplot(X, Y, DX, DY, color='gray')
plt.plot(x_vals,nullclineX,color='red')
plt.plot(x_vals,nullclineY,color='blue')
plt.xlim([-2, 2])
plt.ylim([-2, 2])
plt.title('Phase Plane')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#%% Fig 6.3.1
def xdot(t, x):
    return np.array([-x[0]+x[0]**3,-2*x[1]])

# create grid for phase portrait
x_vals = np.linspace(-2, 2, 100)
y_vals = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

# plotting
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlim([-2, 2])
plt.ylim([-2, 2])
plt.title('Phase Plane')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#%% Fig 6.4.7 Rabbits and Sheep
def xdot(t, x):
    return np.array([x[0]*(3-x[0]-2*x[1]),x[1]*(2-x[0]-x[1])])

# create grid for phase portrait
x_vals = np.linspace(0, 4, 100)
y_vals = np.linspace(0, 4, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

# nullclines solved in wolframalpha- note these don't match the phase portrait
nullclineX1 = (3-x_vals)/2
nullclineY1 = 2-x_vals
nullclineY2 = x_vals*0

# plotting
plt.plot(x_vals,nullclineX1)
plt.plot(x_vals,nullclineY1)
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlim([0, 4])
plt.ylim([0, 4])
plt.title('Phase Plane for Rabbits & Sheep')
plt.xlabel('Rabbits')
plt.ylabel('Sheep')
plt.show()

#%% Fig 6.5.1 Particle in double-well potential
def xdot(t, x):
    return np.array([x[1],x[0]-x[0]**3])

# create grid for phase portrait
x_vals = np.linspace(-2, 2, 100)
y_vals = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

#plotting
plt.streamplot(X, Y, DX, DY, color='gray')
plt.title('Phase Plane for Double-Well Potential')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#%% Fig 6.6.4 nonlinear center- reversible systems
def xdot(t, x):
    return np.array([x[1]-x[1]**3,-x[0]-x[1]**2])

# create grid for phase portrait
x_vals = np.linspace(-2, 2, 100)
y_vals = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

# nullclines solved in wolframalpha- note these don't match the phase portrait
nullclineX1 = -1*np.linspace(1,1,100) # 0 = y-y^3
nullclineX2 = 0*np.linspace(1,1,100) # 0 = y-y^3
nullclineX3 = 1*np.linspace(1,1,100) # 0 = y-y^3
#nullclineY1 = # 0 = -x-y^2

#plotting
plt.streamplot(X, Y, DX, DY, color='gray')
plt.plot(x_vals,nullclineX1,x_vals,nullclineX2,x_vals,nullclineX3)
plt.title('Phase Plane for Nonlinear Center')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#%% Figure 6.7.3 Nonlinear pendulum
def xdot(t, x):
    return np.array([x[1],-np.sin(x[0])])

# create grid for phase portrait
x_vals = np.linspace(-10, 10, 100)
y_vals = np.linspace(-10, 10, 100)
X, Y = np.meshgrid(x_vals, y_vals)
DX, DY = xdot(0, [X, Y]) # vectorize derivative calculations on the grid

#plotting
plt.streamplot(X, Y, DX, DY, color='gray')
plt.title('Phase Plane for Nonlinear Pendulum')
plt.xlabel('Theta')
plt.ylabel('Velocity')
plt.show()
