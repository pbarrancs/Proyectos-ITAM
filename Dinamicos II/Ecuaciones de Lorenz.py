# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 18:30:32 2019

@author: Ana
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

#Parámertos estándar:
# =============================================================================
# rho = 28.0
# sigma = 10.0
# beta = 8.0 / 3.0 
# =============================================================================

#Alternativa 1:
# =============================================================================
# =============================================================================
# rho = 2.0
# sigma = 2.0
# beta = 1.0 
# =============================================================================
# =============================================================================

#Alternativa 2:
# =============================================================================
# rho = 0.5
# sigma = 10.0
# beta = 8.0 / 3.0 
# =============================================================================

#Alternativa 3:
# =============================================================================
# rho = 99.96
# sigma = 10.0
# beta = 8.0 / 3.0 
# =============================================================================

#Alternativa 4:
rho = 14.00
sigma = 10.0
beta = 8.0 / 3.0 

def f(state, t):
  x, y, z = state  # unpack the state vector
  return sigma * (y - x), x * (rho - z) - y, x * y - beta * z  # derivatives

state0 = [0.1, 0.1, 0.1]  #Estado inicial
state1 = [-0.5, -0.1, 0.5]  #Estado inicial
state2 = [-0.1, 0.1, 0.1]  #Estado inicial
state3 = [0.1, -0.1, 0.1]  #Estado inicial
state4 = [0.1, 0.1, -0.1]  #Estado inicial
t = np.arange(0.0, 40.0, 0.01) #Pasos de tiempo

states = odeint(f, state0, t) #Integramos las ecuaciones
# si Rho > 1
X = [0., np.sqrt(beta*(rho-1)), -np.sqrt(beta*(rho-1))]
Y = [0., np.sqrt(beta*(rho-1)), -np.sqrt(beta*(rho-1))]
Z = [0., rho-1, rho-1]
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(states[:,0], states[:,1], states[:,2])
ax.plot(X,Y,Z, 'ro')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()


