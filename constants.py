import math


# first order taylor
E = 1
B = 1
m = 1
q = 1

# error function
analytic_ry_T = 0
w = q * B / m
analytic_rz_T = E / B * 2 * math.pi / w
analytic_T = analytic_ry_T, analytic_rz_T
