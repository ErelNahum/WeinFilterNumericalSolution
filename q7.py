import numpy as np
import constants as c
import math
import matplotlib.pyplot as plt
import q6


def q7(x_intervals, y_intervals):

    energy = np.linspace(c.q6.E0 - c.q6.lE, c.q6.E0 + c.q6.lE, x_intervals)
    radius = np.linspace(-c.q6.R, c.q6.R, y_intervals)

    for y0 in radius:
        for Ei in energy:
            if q6.q6_graph(10000, Ei, y0)[0]:
                plt.plot(y0 / c.q6.R, math.sqrt(Ei / c.q6.E0) - 1, "ro")
    plt.suptitle(r"$(\frac{y_0}{R}, \frac{\delta v}{v_0})$")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    q7(100, 100)
