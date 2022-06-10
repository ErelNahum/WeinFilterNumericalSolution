import constants as c
import numpy as np
import math
import matplotlib.pyplot as plt
import q3


def ay(vz):
    return c.q / c.m * (c.E - c.B * vz)


def az(vy):
    return c.q / c.m * c.B * vy


def midpoint(intervals, draw=False):
    """drawing the graph from t=0 to t=2pi/w"""

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
        f.suptitle("midpoint", fontsize=16)

        axis[0].set_title(r"$(r_y, r_z)$")
        axis[0].grid()
        axis[0].plot(ry, rz)

        axis[1].set_title(r"$(v_y, v_z)$")
        axis[1].grid()
        axis[1].plot(vy, vz)

        plt.show()

    return ry[intervals - 1], rz[intervals - 1]


def runge_kutta(intervals, draw=False):
    """drawing the graph from t=0 to t=T with runge-kutta method"""

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
        k1_vy = ay(vz[t]) * dt
        k1_vz = az(vy[t]) * dt

        k2_vy = ay(vz[t] + 0.5 * k1_vz) * dt
        k2_vz = az(vy[t] + 0.5 * k1_vy) * dt

        k3_vy = ay(vz[t] + 0.5 * k2_vz) * dt
        k3_vz = az(vy[t] + 0.5 * k2_vy) * dt

        k4_vy = ay(vz[t] + k3_vz) * dt
        k4_vz = az(vy[t] + k3_vy) * dt

        vy[t + 1] = vy[t] + (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy) / 6
        vz[t + 1] = vz[t] + (k1_vz + 2 * k2_vz + 2 * k3_vz + k4_vz) / 6

        k1_ry = vy[t] * dt
        k1_rz = vz[t] * dt

        k2_ry = (vy[t] + 0.5 * k1_vy) * dt
        k2_rz = (vz[t] + 0.5 * k1_vz) * dt

        k3_ry = (vy[t] + 0.5 * k2_vy) * dt
        k3_rz = (vz[t] + 0.5 * k2_vz) * dt

        k4_ry = (vy[t] + k3_vy) * dt
        k4_rz = (vz[t] + k3_vz) * dt

        ry[t + 1] = ry[t] + (k1_ry + 2 * k2_ry + 2 * k3_ry + k4_ry) / 6
        rz[t + 1] = rz[t] + (k1_rz + 2 * k2_rz + 2 * k3_rz + k4_rz) / 6

    if draw:
        f, axis = plt.subplots(2, 1)
        f.suptitle("Runge Kutta", fontsize=16)
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
    time = np.zeros(intervals)
    taylor_values = np.zeros(intervals)
    midpoint_values = np.zeros(intervals)
    runge_values = np.zeros(intervals)

    for t in range(intervals - 1):
        time[t + 1] = time[t] + dt
        intervals = math.ceil(c.T / time[t + 1]) + 1
        taylor_values[t + 1] = error(c.analytic_T, q3.q3_graph(intervals))
        midpoint_values[t + 1] = error(c.analytic_T, midpoint(intervals))
        runge_values[t + 1] = error(c.analytic_T, runge_kutta(intervals))

    time = time[1:]
    taylor_values = taylor_values[1:]
    midpoint_values = midpoint_values[1:]
    runge_values = runge_values[1:]

    plt.plot(time, taylor_values, label="taylor")
    plt.plot(time, midpoint_values, label="midpoint")
    plt.plot(time, runge_values, label="runge kutta")
    plt.suptitle("log log - error for dt")
    plt.xscale("log")
    plt.xlabel("dt")
    plt.yscale("log")
    plt.ylabel("error")
    plt.grid()
    plt.legend()

    plt.show()


if __name__ == "__main__":
    midpoint_error_graph(0.0001, 1)


