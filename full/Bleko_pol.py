def forSOE (rho, y1, y2):
    t1 = y2 ** 2
    t2 = -1 + beta
    t6 = alpha - 0.1e1 / 0.2e1
    t10 = k ** 2
    t11 = t10 * rho
    t12 = y1 ** 2
    t14 = rho ** 2
    t15 = 0.1e1 / t14
    t17 = t1 * t12 * t10 * t15
    t19 = t17 ** (2 * alpha)
    t23 = t1 ** 2
    t28 = xi ** 2
    t34 = t14 ** 2
    t44 = t17 ** alpha
    t51 = (alpha + 0.1e1 / 0.2e1) * beta
    return([y2,-y2 * (2 * t1 * t2 * (y2 * alpha * rho - t6 * y1) * t11 * t19 + (-t23 * y1 * beta * t10 * rho + (t10 * (beta * t14 * t28 + 1) * t12 - t34 * t28 * t2) * t1 * y2 - y1 * t10 * rho * t2) * t44 - 2 * (rho * (alpha * beta + 0.1e1 / 0.2e1) * y2 - t51 * y1) * t1 * t11) / t10 / y1 / (t6 * t1 * t2 * t19 + (-beta * t23 / 2 + 0.3e1 / 0.2e1 * beta - 0.3e1 / 0.2e1) * t44 - t51 * t1) * t15 / 2])
def forBC (
  ya1,
  ya2,
  yb1,
  yb2):
    t2 = 2 - 2 * alpha
    t3 = ya2 ** t2
    t5 = 2 * alpha
    t6 = ya1 ** (-t5)
    t7 = k ** (-t5)
    t9 = rho_0 ** 2
    t10 = 0.1e1 / t9
    t11 = t10 ** (-alpha)
    t15 = 2 + 2 * alpha
    t16 = ya2 ** t15
    t17 = -1 + beta
    t19 = ya1 ** t5
    t20 = k ** t5
    t22 = t10 ** alpha
    t25 = ya2 ** 2
    t26 = t25 ** 2
    t30 = yb2 ** t2
    t32 = yb1 ** (-t5)
    t35 = yb2 ** t15
    t37 = yb1 ** t5
    t40 = yb2 ** 2
    t41 = t40 ** 2
    return([-2 * beta * t11 * t3 * t6 * t7 - 2 * t16 * t17 * t19 * t20 * t22 + 2 * beta * t26 + 2 * beta - 2,-2 * beta * t30 * t32 * t7 - 2 * t17 * t20 * t35 * t37 + 2 * beta * t41 + 2 * beta - 2])
import math

def forIntegrand (y1):
    t2 = rho ** 2
    t4 = y1 ** 2
    t6 = k ** 2
    return(2 * xi * ((1 - beta) * t2 + t4 * beta * t6) * math.pi * rho / t6)
