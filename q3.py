import constants as c
import numpy as np
import math
import matplotlib.pyplot as plt


def q3_graph(intervals, draw=False):
    """drawing the graphs from t=0 to t=2pi/w"""

    dt = c.T / (intervals - 1)

    vy = np.zeros(intervals)
    vz = np.zeros(intervals)
    ry = np.zeros(intervals)
    rz = np.zeros(intervals)

    vy[0] = 0
    ry[0] = 0
    vz[0] = 3 * c.E / c.B
    rz[0] = 0

    for t in range(intervals - 1):
        vy[t + 1] = vy[t] + (c.q / c.m * (c.E - c.B * vz[t])) * dt
        vz[t + 1] = vz[t] + c.q * c.B / c.m * vy[t] * dt
        ry[t + 1] = ry[t] + vy[t] * dt
        rz[t + 1] = rz[t] + vz[t] * dt

    if draw:
        f, axis = plt.subplots(2, 1)
        f.suptitle("first order taylor", fontsize=16)

        axis[0].set_title(r"$(r_y, r_z)$")
        axis[0].grid()
        axis[0].plot(ry, rz)

        axis[1].set_title(r"$(v_y, v_z)$")
        axis[1].grid()
        axis[1].plot(vy, vz)

        plt.show()
    return ry[intervals - 1], rz[intervals - 1]


if __name__ == "__main__":
    q3_graph(100000, True)
