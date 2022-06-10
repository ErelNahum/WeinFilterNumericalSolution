import numpy as np
import constants as c
import math
import matplotlib.pyplot as plt


def calculate_v0(E, m):
    E *= c.q6.MEV_TO_J
    return math.sqrt(2 * E / m)


def ay(vz, E):
    return c.q6.q / c.q6.m * (E - c.q6.B * vz)


def az(vy):
    return c.q6.q / c.q6.m * c.q6.B * vy


def q6_graph(intervals, Ei, y0, draw=False):
    collision = False
    v0 = calculate_v0(c.q6.E0, c.q6.m)
    E = v0 * c.q6.B
    v0 = calculate_v0(Ei, c.q6.m)
    dt = c.q6.T / (intervals - 1)

    vy = [0]
    vz = [v0]
    ry = [y0]
    rz = [0]

    for t in range(intervals - 1):
        vy.append(vy[t] + ay(vz[t], E) * dt)
        vz.append(vz[t] + az(vy[t]) * dt)
        ry.append(ry[t] + vy[t] * dt)
        rz.append(rz[t] + vz[t] * dt)

        if abs(ry[t + 1]) > c.q6.R:
            collision = True
            break

        if rz[-1] > 1:
            break

    if draw:
        plt.suptitle(r"$(r_y, r_z)$ : route of a proton with $y_0=-1.5\cdot10^{-4}$ $E_i=4.98\;MeV$")
        plt.grid()
        plt.plot(ry, rz, color="red")
        plt.axvline(x=c.q6.R)
        plt.axvline(x=-c.q6.R)


        # axis[1].set_title(r"$(v_y, v_z)$")
        # axis[1].grid()
        # axis[1].plot(vy, vz)
        plt.show()
    return collision, math.sqrt(vy[-1]**2 + vz[-1]**2)


if __name__ == "__main__":
    print(q6_graph(1000, 5.02, 0, True))


