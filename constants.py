import math
import scipy.constants as s


# first order taylor
E = 1
B = 1
m = 1
q = 1

w = q * B / m
T = 2 * math.pi / w

# error function
analytic_ry_T = 0
analytic_rz_T = E / B * T
analytic_T = analytic_ry_T, analytic_rz_T

class q6:
    m = s.m_p
    q = s.e
    B = 0.5
    w = q * B / m
    T = 2 * math.pi / w
    MEV_TO_J = 1.60217663e-13
    E0 = 5
    R = 3e-3
    lE = 0.25
