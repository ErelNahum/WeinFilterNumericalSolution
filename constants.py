import math


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
