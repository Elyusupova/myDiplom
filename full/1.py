import numpy as np
import Bleko_pol
from scipy.integrate import solve_bvp
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as plt
import math as math


N = 1000
N_plot=100

Bleko_pol.k = 1

rho_0 = 0.1
Bleko_pol.rho_0 = rho_0


xi_s=0.01
xi_f=2.5
xi = 0.3
Bleko_pol.xi = xi


beta = 0.3
Bleko_pol.beta=beta
alpha = 0.25
Bleko_pol.alpha = alpha

rho = np.linspace(rho_0, 1, N)
Bleko_pol.rho=rho

def f(rho,y):
    return np.array (Bleko_pol.forSOE(rho, y[0], y[1]))


def bc (ya, yb):
    return np.array (Bleko_pol.forBC(ya[0], ya[1],yb[0], yb[1]))
                     

def Mom(xi):
    Bleko_pol.xi = xi
    y1= rho
    y2 = np.ones(N)
    y_approx = np.array([y1, y2])
    res_a = solve_bvp(f, bc, rho, y_approx)
    y_plot_a = res_a.sol(rho)[0]
    
    y = list(map(Bleko_pol.forIntegrand, y_plot_a))
    return (integrate.simps(y,rho))

print(Mom(0.3))

import matplotlib.pyplot as plt
points_xi = np.linspace(xi_s,xi_f,N_plot)
points_M = list(map(Mom, points_xi))
plt.plot(points_xi, points_M)
plt.ylabel('Mom')
plt.xlabel('xi ')
plt.show()


#Дима - лох