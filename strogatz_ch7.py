# strogatz_ch7: limit cycles
# Jeff Saucerman 12/30/2025
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

#%% Fig 7.1.1 simple limit cycle, dr/dt = r*(1-r^2); dtheta/dt = 1
def xdot(t, x):
    return np.array([x[0]*(1-x[0]**2),np.ones_like(x[0])])  # note use of ones_like function for vector inputs

rrange = np.linspace(-2,2,100)
R, T = np.meshgrid(rrange, rrange)
DR, DT = xdot(0,[R,T]) 

# plotting phase portrait
plt.streamplot(R, T, DR, DT, color='gray')
plt.xlabel('r')
plt.ylabel('theta')
plt.show()

# Figure 7.1.3 plot one trajectory
x0 = [2, 0.0]  # initial conditions for r, theta
t_span = (0, 10)         
t_eval = np.linspace(0, 10, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval)   

# plotting trajectory in polar coordinates
plt.plot(sol.t, sol.y[0],label='r')
plt.plot(sol.t, sol.y[1],label='theta')
plt.xlabel('t')
plt.legend()
plt.show()

# plotting trajectory in cartesian coordinates
x = sol.y[0]*np.cos(sol.y[1])
y = sol.y[0]*np.sin(sol.y[1])
plt.plot(sol.t,x,label='x')
plt.plot(sol.t,y,label='y')
plt.legend()
plt.show()

#%% Figure 7.1.4 and 7.1.5, van der Pol oscillator
# dx/dt = y; dy/dt = mu*(1-x^2)*y - y
def xdot(t, x, mu):
    return np.array([x[1], mu*(1-x[0]**2)*x[1]-x[0]])

# plot single trajectory
mu = 1.5
x0 = [3, 0.0]  # initial conditions for x, y
t_span = (0, 20)         
t_eval = np.linspace(0, 20, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval, args = (mu,))   
plt.plot(sol.t,sol.y[0],label='x')
plt.plot(sol.t,sol.y[1],label='y')
plt.xlabel('t')
plt.title('Figure 7.1.5 van der Pol oscillator')
plt.legend()
plt.show()

# plot phase portrait with single trajectory
xrange = np.linspace(-4,4,100)
X, Y = np.meshgrid(xrange, xrange)
mu = 1.5
DX, DY = xdot(0,[X,Y], mu) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.plot(sol.y[0],sol.y[1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Figure 7.1.4 van der Pol oscillator')
plt.show()

# Figure 7.5.1 nullclines
nullclineX = 0 # 0 = y, solve for x
nullclineY = 0 # 0 = mu*(1-x^2)*y - y, solve for y
# BUG: nullclineY is supposed to be cubic, confused

#%% Example 7.3.2- glycolytic oscillations
# dx/dt = -x+a*y+x^2*y; dy/dt = b-a*y-x^2*y
def xdot(t, x, a, b):
    return np.array([-x[0]+a*x[1]+x[0]**2*x[1], b-a*x[1]-x[0]**2*x[1]])

# plot single trajectory
a = 0.08
b = 0.6
x0 = [1.0,1.0]  # initial conditions for x, y
t_span = (0, 20)         
t_eval = np.linspace(0, 20, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval, args = (a,b,))   
plt.plot(sol.t,sol.y[0],label='x')
plt.plot(sol.t,sol.y[1],label='y')
plt.xlabel('t')
plt.title('Trajectory from glycolytic oscillator')
plt.legend()
plt.show()

# plot nullclines
xrange = np.linspace(0,5,100)
nullclineX = xrange/(a+xrange**2) # y = x/(a+x^2)
nullclineY = b/(a+xrange**2)     # y = b/(a+x^2)
plt.plot(xrange,nullclineX,label='nullclineX')
plt.plot(xrange,nullclineY,label='nullclineY')
plt.title('Figure 7.3.4: nullclines')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

# Figure 7.3.8, phase portrait for a = 0.08, b = 0.6
X, Y = np.meshgrid(xrange, xrange)
DX, DY = xdot(0,[X,Y], a, b) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.plot(sol.y[0],sol.y[1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Figure 7.3.8 glycolytic oscillator')
plt.show()
