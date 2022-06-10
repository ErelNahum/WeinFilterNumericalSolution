import q6
import random
import constants as c
import matplotlib.pyplot as plt


def q8(n):
    result = []
    for _ in range(n):
        if _ % 1000 == 0:
            print(_)
        Ei = random.random() * 2 * c.q6.lE + c.q6.E0 - c.q6.lE
        y0 = random.random() * 2 * c.q6.R - c.q6.R
        x = q6.q6_graph(100000, Ei, y0)
        if not x[0]:
            result.append(x[1])

    plt.hist(result, bins=70)
    plt.xlabel(r"exit velocity[$\frac{m}{s}$]")
    plt.suptitle("exit velocities distribution")
    plt.show()


if __name__ == '__main__':
    q8(100000)
