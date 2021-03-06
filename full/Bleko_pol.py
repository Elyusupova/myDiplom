def forSOE (rho, y1, y2):
    t3 = alpha - 0.1e1 / 0.2e1
    t7 = -1 + beta
    t9 = k ** 2
    t10 = y2 ** 2
    t11 = t9 * t10
    t12 = y1 ** 2
    t14 = rho ** 2
    t15 = 0.1e1 / t14
    t17 = t10 * t12 * t9 * t15
    t19 = t17 ** (2 * alpha)
    t23 = t10 ** 2
    t28 = xi ** 2
    t34 = t14 ** 2
    t44 = t17 ** alpha
    t50 = alpha + 0.1e1 / 0.2e1
    return([y2,-(2 * rho * (y2 * alpha * rho - t3 * y1) * t7 * t11 * t19 + (-t23 * y1 * beta * t9 * rho + (t9 * (beta * t14 * t28 + 1) * t12 - t34 * t28 * t7) * t10 * y2 - y1 * t9 * rho * t7) * t44 - 2 * rho * (rho * (alpha * beta + 0.1e1 / 0.2e1) * y2 - y1 * t50 * beta) * t11) * y2 * t15 / (t7 * t10 * t3 * t19 + (-beta * t23 / 2 + 0.3e1 / 0.2e1 * beta - 0.3e1 / 0.2e1) * t44 - beta * t10 * t50) / t9 / y1 / 2])
def forBC (
  ya1,
  ya2,
  yb1,
  yb2):
    return([-2 * beta * ya2 ** (2 - 2 * alpha) * ya1 ** (-2 * alpha) * k ** (-2 * alpha) * (1 / rho_0 ** 2) ** (-alpha) - 2 * ya2 ** (2 + 2 * alpha) * (-1 + beta) * ya1 ** (2 * alpha) * k ** (2 * alpha) * (1 / rho_0 ** 2) ** alpha + 2 * beta * ya2 ** 4 + 2 * beta - 2,-2 * beta * yb2 ** (2 - 2 * alpha) * yb1 ** (-2 * alpha) * k ** (-2 * alpha) - 2 * yb2 ** (2 + 2 * alpha) * (-1 + beta) * yb1 ** (2 * alpha) * k ** (2 * alpha) + 2 * beta * yb2 ** 4 + 2 * beta - 2])
import math

def forIntegrand (y1):
    t2 = rho ** 2
    t4 = y1 ** 2
    t6 = k ** 2
    return(2 * xi * ((1 - beta) * t2 + t4 * beta * t6) * math.pi * rho / t6)
