# strogatz_ch8.py: Bifurcations revisited
# Jeff Saucerma, 12/30/2025
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

#%% Fig 8.1.1 saddlenode bifurcation
# dx/dt = mu-x^2; dy/dt = -y
def xdot(t, x, mu):
    return np.array([mu-x[0]**2,-x[1]]) 
 
xrange = np.linspace(-4,4,100)
X, Y = np.meshgrid(xrange, xrange)
mu = 1
DX, DY = xdot(0,[X,Y], mu) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'Fig 8.1.1 saddle node $\mu$ = {mu}')

# Figure 8.1.2 nullclines
nullclineX = mu*np.sqrt(xrange) # 0 = mu-x^2, solve for y?
nullclineY = 0 # 0 = -y, solve for y
# BUG? nullclines don't match the text although it notes the sketch
# is for a general system, not the same as before

#%% Example 8.1.1 genetic toggle switch Griffith (1971)
# dx/dt = -a*x+y; dy/dt = x^2/(1+x^2) - b*y
def xdot(t, x, a, b):
    return np.array([-a*x[0]+x[1], x[0]**2/(1+x[0]**2) - b*x[1]]) 
a = 1
b = 0.4
# plot nullclines
xrange = np.linspace(0,2,100)
nullclineX = a*xrange                     # 0 = -a*x+y, solve for y
nullclineY = xrange**2/(1+xrange**2)/b    # 0 = x^2/(1+x^2) - b*y, solve for y
plt.plot(xrange,nullclineX,label='nullclineX')
plt.plot(xrange,nullclineY,label='nullclineY')
plt.title('Figure 8.1.4: nullclines')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

# plot phase portrait
X, Y = np.meshgrid(xrange, xrange)
DX, DY = xdot(0,[X,Y], a, b) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'Fig 8.1.5 saddle node b = {b}')

#%% Example 8.1.2 supercritical pitchfork bifurcation
# dx/dt = mu*x-x^3; dy/dt = -y
def xdot(t, x, mu):
    return np.array([mu*x[0] - x[0]**3, -x[1]]) 

# plot phase portrait
xrange = np.linspace(-2,2,100)
X, Y = np.meshgrid(xrange, xrange)
mu = 1
DX, DY = xdot(0,[X,Y], mu) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'Fig 8.1.6 supercritical pitchfork $\mu$ = {mu}')

#%% transcritical bifurcation
# dx/dt = mu*x-x^2; dy/dt = -y
def xdot(t, x, mu):
    return np.array([mu*x[0] - x[0]**2, -x[1]]) 

# plot phase portrait
xrange = np.linspace(-2,2,100)
X, Y = np.meshgrid(xrange, xrange)
mu = 1
DX, DY = xdot(0,[X,Y], mu) 
plt.streamplot(X, Y, DX, DY, color='gray')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'transcritical bifurcation $\mu$ = {mu}')

#%% supercritical Hopf bifurcation
# dr/dt = mu*r-r^3; dtheta/dt = omega+b*r^2
def xdot(t, x, mu, omega, b):
    return np.array([mu*x[0]-x[0]**3, omega + b*x[0]**2])

# plot single trajectory
mu = -1
omega = 10
b = 1
x0 = [0.5,0.5]  # initial conditions for x, y
t_span = (0, 20)         
t_eval = np.linspace(0, 20, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval, args = (mu,omega,b,))   
x = sol.y[0]*np.cos(sol.y[1])
y = sol.y[0]*np.sin(sol.y[1])
plt.plot(x,y)
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Phase plane for $\mu$ = {mu}')
plt.legend()
plt.show()

#%% subcritical Hopf bifurcation
# dr/dt = mu*r+r^3-r^5; dtheta/dt = omega+b*r^2
def xdot(t, x, mu, omega, b):
    return np.array([mu*x[0]+x[0]**3-x[0]**5, omega + b*x[0]**2])

# plot single trajectory
mu = 1
omega = 10
b = 1
x0 = [0.5,0.5]  # initial conditions for x, y
t_span = (0, 20)         
t_eval = np.linspace(0, 20, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval, args = (mu,omega,b,))   
x = sol.y[0]*np.cos(sol.y[1])
y = sol.y[0]*np.sin(sol.y[1])
plt.plot(x,y)
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Phase plane for $\mu$ = {mu}')
plt.show()

#%% Belousov-Zhabotinsky oscillating reaction
# dx/dt = a-x-4*x*y/(1+x^2); dy/dt = b*x*(1-y/(1+x^2))
def xdot(t, x, a, b):
    return np.array([a-x[0]-4*x[0]*x[1]/(1+x[0]**2), b*x[0]*(1-x[1]/(1+x[0]**2))])

# plot single trajectory
a = 10
b = 4 # stable fixed point for b = 4, limit cycle for b = 2
x0 = [0.5,0.5]  # initial conditions for x, y
t_span = (0, 20)         
t_eval = np.linspace(0, 20, 500)
sol = scipy.integrate.solve_ivp(xdot, t_span, x0, t_eval=t_eval, args = (a,b,))   
plt.plot(sol.t,sol.y[0],label='x')
plt.plot(sol.t,sol.y[1],label='y')
plt.xlabel('t')
plt.legend()
plt.title(f'BZ reaction, a = {a}, b = {b}')
plt.show()

plt.plot(sol.y[0],sol.y[1])
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'BZ reaction, a = {a}, b = {b}')
plt.show()