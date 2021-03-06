import numpy as np
import Bleko_upr
from scipy.integrate import solve_bvp
import scipy.integrate as integrate
import matplotlib as plt
from scipy.integrate import quad
import math as math




N = 1000
N_plot=100

Bleko_upr.k = 1

rho_0 = 0.9
Bleko_upr.rho_0 = rho_0

Bleko_upr.pi=math.pi
xi_s=0.01
xi_f=2.5
xi = 0.3
Bleko_upr.xi = xi

alpha = 0.15
Bleko_upr.alpha = alpha

rho = np.linspace(rho_0, 1, N)
Bleko_upr.rho=rho

def f(rho,y):
    return np.array (Bleko_upr.forSOE(rho, y[0], y[1]))


def bc (ya, yb):
    return np.array (Bleko_upr.forBC(ya[0], ya[1],yb[0], yb[1]))
                     


def Mom(xi):
    Bleko_upr.xi = xi
    y1= rho
    y2 = np.ones(N)
    y_approx = np.array([y1, y2])
    res_a = solve_bvp(f, bc, rho, y_approx)
    y_plot_a = res_a.sol(rho)[0]

    y = list(map(Bleko_upr.forIntegrand, y_plot_a))
    return (integrate.simps(y,rho))

print(Mom(0.3))

import matplotlib.pyplot as plt
points_xi = np.linspace(xi_s,xi_f,N_plot)
points_M = list(map(Mom, points_xi))
plt.plot(points_xi, points_M)
plt.ylabel('Mom')
plt.xlabel('xi ')
plt.show()


