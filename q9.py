import random
import constants as c
import q6
import numpy as np

def get_pass_rate(n):
    result = 0
    for _ in range(n):
        if _ % 10000 == 0:
            print(f"progress: {_}")
        Ei = random.random() * 2 * c.q6.lE + c.q6.E0 - c.q6.lE
        y0 = random.random() * 2 * c.q6.R - c.q6.R
        if not q6.q6_graph(1000, Ei, y0)[0]:
            result += 1
    return result / n


if __name__ == '__main__':
    result = []
    for _ in range(100):
        result.append(get_pass_rate(10000))
    print(np.average(result))
    print(np.std(result))

