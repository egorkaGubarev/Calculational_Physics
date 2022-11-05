import numpy as np
import matplotlib.pyplot as plt


def get_level(path, level_number_l):
    with open(path, 'r') as file:
        lines = file.readlines()
        level = np.array(lines[level_number_l].split(), dtype=float)
    return level


def get_scaling_decomposition(coefficients_l, points_l, level):
    amount = len(coefficients_l)
    points_per_gap = points_l / amount
    y_l = np.zeros(points_l)

    for point_number in range(points_l):
        gap_number = int(point_number / points_per_gap)
        y_l[point_number] = coefficients_l[gap_number] * 2 ** ((level - 10) / 2)

    return y_l


path_scaling = '../../../VSProjects/Wavelet/Wavelet/scaling.txt'
path_wavelet = '../../../VSProjects/Wavelet/Wavelet/wavelet.txt'

left_x = 0
right_x = 1

levels = 10
points = 1000

x = np.linspace(left_x, right_x, points)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
plt.tight_layout()

coefficients = get_level(path_scaling, 0)
y = get_scaling_decomposition(coefficients, points, 0)
ax1.set_title('0 level')
ax1.set_xlabel('x')
ax1.set_ylabel('f')
ax1.plot(x, y)

coefficients = get_level(path_scaling, 2)
y = get_scaling_decomposition(coefficients, points, 2)
ax2.set_title('2 level')
ax2.set_xlabel('x')
ax2.set_ylabel('f')
ax2.plot(x, y)

coefficients = get_level(path_scaling, 4)
y = get_scaling_decomposition(coefficients, points, 4)
ax3.set_title('4 level')
ax3.set_xlabel('x')
ax3.set_ylabel('f')
ax3.plot(x, y)

coefficients = get_level(path_scaling, 10)
y = get_scaling_decomposition(coefficients, points, 10)
ax4.set_title('10 level')
ax4.set_xlabel('x')
ax4.set_ylabel('f')
ax4.plot(x, y)


plt.show()
