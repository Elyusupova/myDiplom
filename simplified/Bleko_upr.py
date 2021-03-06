def forSOE (rho, y1, y2):
    t1 = y2 ** 2
    t3 = rho ** 2
    t6 = y1 ** 2
    t7 = t6 * y1
    return([y2,-y2 * (rho * t1 * t3 * y2 - t7) / t7 / rho / 3])
def forBC (
  ya1,
  ya2,
  yb1,
  yb2):
    t1 = ya2 ** 2
    t7 = yb2 ** 2
    return([2 * k * t1 * ya1 * ya2 - 2 * rho_0,2 * k * t7 * yb1 * yb2 - 2])
import math

def forIntegrand (y1):
    t1 = rho ** 2
    t4 = k ** 2
    return(4 * math.pi * t1 * rho * xi / t4)
