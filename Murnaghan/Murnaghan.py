def forSOE (rho, y1, y2):
    t1 = rho ** 2
    t2 = t1 ** 2
    t4 = 2 * m
    t5 = l + t4
    t7 = y2 ** 2
    t8 = t7 ** 2
    t13 = xi ** 2
    t14 = t13 * t1
    t15 = t14 + 1
    t22 = y1 ** 2
    t24 = k ** 2
    t25 = t24 * l
    t26 = 3 * l
    t27 = 2 * mu
    t29 = t1 * (t25 - t26 + lam - t4 + t27)
    t37 = t15 ** 2
    t38 = l * t37
    t40 = t24 - 3
    t41 = t40 * l
    t42 = n / 2
    t45 = t13 * (t41 + m - t42 + lam) * t1
    t46 = -m + t42
    t49 = t1 * (t24 * t46 + lam + m + t41 - t42 + t45)
    t54 = l * t15
    t56 = t22 ** 2
    t68 = t24 ** 2
    t70 = (-6 * t24 + t68 + 0.9e1 / 0.4e1) * l
    t71 = 2 * lam
    t73 = (t71 + t4 - n) * t24
    t74 = n / 4
    t75 = 0.3e1 / 0.2e1 * lam
    t77 = t2 * (t70 + t73 - mu + t74 - t75)
    return([y2,(-4 * t2 * rho * t5 * t8 * y2 - 12 * y1 * l * t2 * t15 * t8 - (8 * (t14 * l - l) * t22 + 8 * t29) * t1 * rho * t7 * y2 - 8 * t1 * y1 * (t22 * t38 + t49) * t7 - 4 * rho * (t54 * (t14 - 3) * t56 + 2 * t1 * (-t40 * l - t24 * t46 - lam - m + t42 + t45) * t22 + t77) * y2 + 4 * y1 * (t37 * t15 * t5 * t56 + 2 * t15 * t1 * ((2 * m * t24 + lam + t27 - t4 + t41) * t13 * t1 + t25 - t26 - t4 + t27 + lam) * t22 + t2 * ((t70 + 2 * t68 * m + (t71 - 4 * m + 4 * mu) * t24 - mu + t74 - t75) * t13 * t1 + t70 + t73 - mu + t74 - t75))) / t1 / (5 * t2 * t5 * t8 + 6 * t1 * (t22 * t54 + t29) * t7 + t38 * t56 + 2 * t49 * t22 + t77) / 4])
def forBC (
  ya1,
  ya2,
  yb1,
  yb2):
    t1 = xi ** 2
    t2 = rho_0 ** 2
    t4 = t1 * t2 + 1
    t5 = t4 ** 2
    t7 = ya1 ** 2
    t8 = t7 ** 2
    t11 = ya2 ** 2
    t14 = k ** 2
    t15 = t14 - 3
    t16 = l * t15
    t17 = n / 2
    t26 = 2 * m
    t27 = 2 * lam
    t31 = l + t26
    t32 = t11 ** 2
    t39 = 2 * t14 * l - 6 * l - 4 * m + 4 * mu + t27
    t42 = t14 ** 2
    t44 = (-6 * t14 + t42 + 0.9e1 / 0.4e1) * l
    t45 = t27 + t26 - n
    t46 = t45 * t14
    t47 = n / 4
    t48 = 0.3e1 / 0.2e1 * lam
    t50 = t2 ** 2
    t55 = t1 + 1
    t56 = t55 ** 2
    t58 = yb1 ** 2
    t59 = t58 ** 2
    t61 = 2 * t55
    t63 = yb2 ** 2
    t72 = t63 ** 2
    return([(4 * l * t5 * t8 + 4 * (2 * l * t4 * t11 + 2 * t1 * (t16 + m - t17 + lam) * t2 + 2 * t16 + 2 * (-m + t17) * t14 + t26 - n + t27) * t2 * t7 + 4 * (t11 * t39 + t31 * t32 - mu + t44 + t46 + t47 - t48) * t50) * ya2,4 * yb2 * (l * t56 * t59 + (t61 * l * t63 + t61 * t15 * l + (-t26 + n) * t14 + t45 * t55) * t58 + t31 * t72 + t39 * t63 + t44 + t46 + t47 - mu - t48)])
import math

def forIntegrand (y1, y2):
    t1 = y1 ** 2
    t4 = rho ** 2
    t5 = t4 ** 2
    t6 = y2 ** 2
    t7 = t6 ** 2
    t10 = xi ** 2
    t13 = k ** 2
    t14 = t13 * l
    t24 = t10 ** 2
    t25 = 2 * m
    t26 = l + t25
    t28 = t1 ** 2
    t33 = 2 * t13 * m
    t40 = t13 ** 2
    t45 = 2 * lam
    t46 = 4 * m
    t47 = 4 * mu
    return([2 * xi * t1 * math.pi * (t5 * t7 * l + (2 * (l * t10 * t1 + t14 - 3 * l + m - n / 2 + lam) * t4 + 2 * l * t1) * t4 * t6 + (t24 * t26 * t28 + 2 * t10 * ((t13 - 3) * l + t33 - t25 + 2 * mu + lam) * t1 + (-6 * t13 + t40 + 0.9e1 / 0.4e1) * l + 2 * t40 * m + (t45 - t46 + t47) * t13 - mu + n / 4 - 0.3e1 / 0.2e1 * lam) * t5 + (2 * t1 * t10 * t26 - 6 * l + 2 * t14 + t33 + t45 - t46 + t47) * t1 * t4 + t28 * t26) / t4 / rho / mu])
