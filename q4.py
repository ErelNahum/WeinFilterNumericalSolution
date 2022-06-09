import constants as c
import numpy as np
import math
import matplotlib.pyplot as plt
import q3


def ay(vz):
    return c.q / c.m * (c.E - c.B * vz)


def az(vy):
    return c.q / c.m * c.B * vy


def q4_midpoint(dt, draw=False):
    """drawing the graph from t=0 to t=2pi/w"""
    w = c.q * c.B / c.m
    T = 2 * np.pi / w

    intervals = math.ceil(T / dt)

    vy = np.zeros(intervals)
    vz = np.zeros(intervals)
    ry = np.zeros(intervals)
    rz = np.zeros(intervals)

    vy[0] = 0
    ry[0] = 0
    vz[0] = 3 * c.E / c.B
    rz[0] = 0

    for t in range(intervals - 1):
        vy_k1 = dt * az(vy[t])
        vy_k2 = dt * ay(vz[t] + 0.5 * vy_k1)
        vy[t + 1] = vy[t] + vy_k2

        vz_k1 = dt * ay(vz[t])
        vz_k2 = dt * az(vy[t] + 0.5 * vz_k1)
        vz[t + 1] = vz[t] + vz_k2

        ry_k1 = dt * ay(vz[t])
        ry_k2 = dt * (vy[t] + 0.5 * ry_k1)
        ry[t + 1] = ry[t] + ry_k2

        rz_k1 = dt * az(vy[t])
        rz_k2 = dt * (vz[t] + 0.5 * rz_k1)
        rz[t + 1] = rz[t] + rz_k2

    if draw:
        f, axis = plt.subplots(2, 1)

        axis[0].set_title(r"$(r_y, r_z)$")
        axis[0].grid()
        axis[0].plot(ry, rz)

        axis[1].set_title(r"$(v_y, v_z)$")
        axis[1].grid()
        axis[1].plot(vy, vz)

        plt.show()

    return ry[intervals - 1], rz[intervals - 1]


def error(analytic, calculation):
    return (analytic[0] - calculation[0]) ** 2 + (analytic[1] - calculation[1]) ** 2


def midpoint_error_graph(dt, max_dt):
    intervals = math.ceil(max_dt / dt)
    x = np.zeros(intervals)
    y = np.zeros(intervals)
    q = np.zeros(intervals)

    for t in range(intervals - 1):
        x[t + 1] = x[t] + dt
        y[t + 1] = error(c.analytic_T, q4_midpoint(x[t+1]))
        q[t + 1] = error(c.analytic_T, q3.q3_graph(x[t+1]))
    x = x[1:]
    y = y[1:]
    q = q[1:]
    plt.plot(x, y, label="midpoint")
    plt.plot(x, q, label="taylor")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show()

q4_midpoint(0.0001, True)
midpoint_error_graph(0.0001, 0.01)