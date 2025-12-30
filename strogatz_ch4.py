# strogatz_ch4: Flows on the circle
import numpy as np
import matplotlib.pyplot as plt
#import scipy.integrate
#import scipy.optimize
#import sympy as sp


#%% Fig 4.3.2. xdot = r + x**2
def thetadot(t,theta,omega,a): 
    return omega - a*np.sin(theta)

omega = 1
arange = np.array([0.5, 1, 1.5])
thetarange = np.arange(-5, 5, 0.1)

fig, ax432 = plt.subplots()
for a in arange:
    ax432.plot(thetarange,thetadot(0,thetarange,omega,a))
plt.xlabel('theta')
plt.ylabel('thetadot')
plt.legend(['a=0.5','a=1','a=1.5'])
plt.show()    
#    sol = scipy.integrate.solve_ivp(thetadot,tspan,[y0],rtol=1e-6)
#    ax213.plot(sol.t,sol.y.T)

