# strogatz_ch3: Bifurcations
# Jeff Saucerman 12/30/2025
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.optimize
import sympy as sp

#%% Fig 3.1.1. xdot = r + x**2
def xdot(t,x,r): 
    return r + x**2
def equation(x,r): 
    return xdot(0,x,r)           # make xdot autonomous for root-finding

r = 1 # -1, 0, 1
sol1 = scipy.optimize.root(equation,-1.1,args=r)    # find two roots based on guesses
sol2 = scipy.optimize.root(equation,1.1,args=r)
# NOTE: sol1.succees = False, indicating no root is found when when r = 1 (equation = 1, not 0)?

xrange = np.linspace(-2,2,100)
fig, ax311 = plt.subplots()
ax311.plot(xrange,xdot(0,xrange,r))
# plot roots (need to know stability by viewing phase plot)
#ax311.plot(sol1.x[0],0,marker='o',markersize=12,markerfacecolor='red') # for r = -1, stable f.p. at x = -1
#ax311.plot(sol2.x[0],0,marker='o',markersize=12,markerfacecolor='white') # for r = -1, unstable f.p. at x = 1
ax311.set_xlabel('x')
ax311.set_ylabel('xdot')
ax311.set_title('Figure 3.1.1')
plt.show()

#%% Fig 3.1.4 bifurcation diagram for xdot = r + x**2
fig, ax314 = plt.subplots()
rrange = np.linspace(-1,1,20)
for r in rrange:
    sol1 = scipy.optimize.root(equation,-2,args=r)    # find two roots based on guesses
    sol2 = scipy.optimize.root(equation,2,args=r)

    if sol1.success:
        ax314.plot(r,sol1.x[0],marker='o',markersize=12,markerfacecolor='red') # For stability/existence see below
    if sol2.success:
        ax314.plot(r,sol2.x[0],marker='o',markersize=12,markerfacecolor='white') # 
ax314.set_xlabel('r')
ax314.set_ylabel('x')
ax314.set_xlim(-1,1)
ax314.set_title('Figure 3.1.4')
plt.show()

#%% Symbolic stability analysis of xdot = r + x**2
x, r = sp.symbols('x r')
f = r + x**2
fixedpoints = sp.solve(f,x)

# stability of fixed points if f'(x*) < 0
fprime = sp.diff(f,x)
fprime_xstar = np.array([fprime.subs(x, y) for y in fixedpoints])
# -2*sqrt(-r) < 0 for r < 0; 2*sqrt(-r) > 1 for r < 0

# Evaluate the fixed points at each r value. (<0 is stable)
# NOTE: the values look right, but I can't extract real parts to plot the vals
stable0 = np.array([fprime_xstar[0].subs(r,ri) for ri in rrange])
stable1 = np.array([fprime_xstar[1].subs(r,ri) for ri in rrange])
# isStable0 = np.real(stable0) < 0 # should work but getting "invalid comparison" error

#%% Fig 3.2.1. xdot = r*x - x**2
def xdot(t,x,r): 
    return r*x - x**2
def equation(x,r): 
    return xdot(0,x,r)           # make xdot autonomous for root-finding

r = 1 # -1, 0, 1
sol1 = scipy.optimize.root(equation,-1.1,args=r)    # find two roots based on guesses
sol2 = scipy.optimize.root(equation,1.1,args=r)
# NOTE: sol1.succees = False, indicating no root is found when when r = 1 (equation = 1, not 0)?

xrange = np.linspace(-1.5,1.5,100)
fig, ax321 = plt.subplots()
ax321.plot(xrange,xdot(0,xrange,r))
# plot roots (need to know stability by viewing phase plot)
ax321.plot(sol1.x[0],0,marker='o',markersize=12,markerfacecolor='white') # for r = -1, unstable f.p. at x = -1
ax321.plot(sol2.x[0],0,marker='o',markersize=12,markerfacecolor='red') # for r = -1, stable f.p. at x = 0
ax321.set_xlabel('x')
ax321.set_ylabel('xdot')
ax321.set_title('Figure 3.2.1 r = 1')
plt.show()

#%% Fig 3.2.2 bifurcation diagram for xdot = r*x - x**2, transcritical bifurcation
fig, ax322 = plt.subplots()
rrange = np.linspace(-1,1,20)
for r in rrange:
    sol1 = scipy.optimize.root(equation,-2,args=r)    # find two roots based on guesses
    sol2 = scipy.optimize.root(equation,2,args=r)
    if sol1.success:
        ax322.plot(r,sol1.x[0],marker='o',markersize=12,markerfacecolor='white') # For stability/existence see below
    if sol2.success:
        ax322.plot(r,sol2.x[0],marker='o',markersize=12,markerfacecolor='red') # 
ax322.set_xlabel('r')
ax322.set_ylabel('x')
ax322.set_xlim(-1,1)
ax322.set_title('Figure 3.2.2')
plt.show()

#%% Fig 3.4.1 super-critical pitchfork bifurcation with xdot = r*x - x^3
def xdot(t,x,r): 
    return r*x - x**3
def equation(x,r): 
    return xdot(0,x,r)           # make xdot autonomous for root-finding
r = 1 # -1, 0, 1
sol1 = scipy.optimize.root(equation,-1.1,args=r)    # find two roots based on guesses
sol2 = scipy.optimize.root(equation,0.1,args=r)
sol3 = scipy.optimize.root(equation,1.1,args=r)
xrange = np.linspace(-1.5,1.5,100)
fig, ax341 = plt.subplots()
ax341.plot(xrange,xdot(0,xrange,r))
# plot roots (need to know stability by viewing phase plot)
ax341.plot(sol1.x[0],0,marker='o',markersize=12,markerfacecolor='red') # for r = -1, unstable f.p. at x = -1
ax341.plot(sol2.x[0],0,marker='o',markersize=12,markerfacecolor='white') # for r = -1, stable f.p. at x = 0
ax341.plot(sol3.x[0],0,marker='o',markersize=12,markerfacecolor='red') # for r = -1, stable f.p. at x = 0
ax341.set_xlabel('x')
ax341.set_ylabel('xdot')
ax341.set_title('Figure 3.4.1 r = 1')
plt.show()

#%% Fig 3.4.2 bifurcation diagram for xdot = r*x - x^3, super-critical pitchfork bifurcation
fig, ax342 = plt.subplots()
rrange = np.linspace(-1,1,20)
for r in rrange:
    sol1 = scipy.optimize.root(equation,-2,args=r)    # find two roots based on guesses
    sol2 = scipy.optimize.root(equation,0,args=r)
    sol3 = scipy.optimize.root(equation,2,args=r)
    if sol1.success:
        ax342.plot(r,sol1.x[0],marker='o',markersize=12,markerfacecolor='red') # For stability/existence see below
    if sol2.success:
        ax342.plot(r,sol2.x[0],marker='o',markersize=12,markerfacecolor='white') # 
    if sol2.success:
        ax342.plot(r,sol3.x[0],marker='o',markersize=12,markerfacecolor='red') # 
ax342.set_xlabel('r')
ax342.set_ylabel('x')
ax342.set_xlim(-1,1)
ax342.set_title('Figure 3.4.2')
plt.show()

#%% Fig 3.4.5 Potential plots for xdot = r*x - x^3
def V(x,r):     # V = -integral(f(x),dx)
    return -0.5*r*x**2 + 0.25*x**4 # V(x) = -integral(xdot,dx)
r = 1 # -1, 0, 1
x = np.linspace(-2,2,100)
fig, ax345 = plt.subplots()
ax345.plot(x,V(x,r))
ax345.set_xlabel('x')
ax345.set_ylabel('V(x,r)')
ax345.set_title('Figure 3.4.5')
plt.show()

#%% Fig 3.4.6 subcritical pitchfork bifurcation diagram with xdot = r*x + x^3
def xdot(t,x,r): 
    return r*x + x**3
def equation(x,r): 
    return xdot(0,x,r)           # make xdot autonomous for root-finding
fig, ax346 = plt.subplots()
rrange = np.linspace(-1,1,20) 
for r in rrange:
    sol1 = scipy.optimize.root(equation,-2,args=r)    # find two roots based on guesses
    sol2 = scipy.optimize.root(equation,0,args=r)
    sol3 = scipy.optimize.root(equation,2,args=r)
    if sol1.success:
        ax346.plot(r,sol1.x[0],marker='o',markersize=12,markerfacecolor='white') # For stability/existence see below
    if sol2.success:
        ax346.plot(r,sol2.x[0],marker='o',markersize=12,markerfacecolor='red') # 
    if sol2.success:
        ax346.plot(r,sol3.x[0],marker='o',markersize=12,markerfacecolor='white') # 
ax346.set_xlabel('r')
ax346.set_ylabel('x')
ax346.set_title('Figure 3.4.6')
plt.show()

#%% 3.7.4 Insect Outbreak
def xdot(t,x,r,k): # default parameter r = 0.5, k = 20
    return r*x*(1-x/k) - x**2/(1 + x**2)
def equation(x,k):
    r = 0.5
    return xdot(0,x,r,k)           # make xdot autonomous for root-finding
xrange = np.linspace(-0.1,2,100)
fig, ax374 = plt.subplots()
ax374.plot(xrange,xdot(0,xrange,r=0.5,k=20))
ax374.set_xlabel('x')
ax374.set_ylabel('xdot')

# plot roots at 1 k value (need to know stability by viewing phase plot)
k = 20
sol1 = scipy.optimize.root(equation,-0.1,args=k)    # find two roots based on guesses
sol2 = scipy.optimize.root(equation,0.5,args=k)
sol3 = scipy.optimize.root(equation,2,args=k)
ax374.plot(sol1.x[0],0,marker='o',markersize=12,markerfacecolor='white') # for k = 20, unstable f.p. at x = -1
ax374.plot(sol2.x[0],0,marker='o',markersize=12,markerfacecolor='red') # for k = 20, stable f.p. at x = 0
ax374.plot(sol3.x[0],0,marker='o',markersize=12,markerfacecolor='white') # for k = 20, stable f.p. at x = 0
ax374.set_title('Figure 3.7.4')
plt.show()

# bifurcation diagram wtih k as bifurcation parameter
fig, ax375 = plt.subplots()
krange = np.linspace(1,40,40) 
for k in krange:
    sol1 = scipy.optimize.root(equation,-0.1,args=k)    # find three roots based on guesses
    sol2 = scipy.optimize.root(equation,0.5,args=k)
    sol3 = scipy.optimize.root(equation,2,args=k)
    print(f"k = {k}, sol1 = {sol1.x[0]}, sol2 = {sol2.x[0]}, sol3 = {sol3.x[0]}")
    if sol1.success:
        ax375.plot(k, sol1.x[0],marker='o',markersize=12,markerfacecolor='white') # for k = 20, unstable f.p. at x = -1
    if sol2.success:
        ax375.plot(k, sol2.x[0],marker='o',markersize=12,markerfacecolor='red') # for k = 20, stable f.p. at x = 0
    if sol3.success:
        ax375.plot(k, sol3.x[0],marker='o',markersize=12,markerfacecolor='white') # for k = 20, stable f.p. at x = 0
ax375.set_xlabel('k')
ax375.set_ylabel('x')
ax375.set_title('Figure 3.7.5')
plt.show()




