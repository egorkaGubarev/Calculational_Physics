import numpy as np
import matplotlib.pyplot as plt

path_init = '../../../VSProjects/Calculational_physics/Trolley/init.txt'
path_res = '../../../VSProjects/Calculational_physics/Trolley/result.txt'


def get_delta(vector):
    size = len(vector)
    size_delta = size - 1
    delta_vector = np.zeros(size_delta)
    for i in range(size_delta):
        next_index = i + 1
        delta = vector[next_index] - vector[i]
        delta_vector[i] = delta
    return delta_vector


init = np.loadtxt(path_init)
res = np.loadtxt(path_res)

delta_init = get_delta(init)
delta_res = get_delta(res)

fig, ((ax_1, ax_2), (ax_3, ax_4), (ax_5, ax_6)) = plt.subplots(3, 2)
plt.tight_layout()
ax_1.hist(delta_init)
ax_1.set_title('Initial')
ax_1.set_xlabel('Distance, m')
ax_1.set_ylabel('Amount')
ax_2.hist(delta_res)
ax_2.set_title('Result')
ax_2.set_xlabel('Distance, m')
ax_2.set_ylabel('Amount')
ax_3.hist(init)
ax_3.set_xlabel('x, m')
ax_3.set_ylabel('Amount')
ax_4.hist(res)
ax_4.set_xlabel('x, m')
ax_4.set_ylabel('Amount')
ax_5.plot(init, np.arange(0, len(init)))
ax_5.set_ylabel('Number of trolley')
ax_5.set_xlabel('x, m')
ax_6.plot(res, np.arange(0, len(res)))
ax_6.set_ylabel('Number of trolley')
ax_6.set_xlabel('x, m')
plt.show()
