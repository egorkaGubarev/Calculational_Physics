import numpy as np
import matplotlib.pyplot as plt

path_res = '../../../VSProjects/Calculational_physics/Trolley_phase_diagram/collapse.txt'

result = np.loadtxt(path_res)
pass_per = result[:, 0]
time = result[:, 1]
plt.plot(pass_per, time, '.')
plt.title('Phase Diagram')
plt.xlabel('Period of attempt to generate passenger, s')
plt.ylabel('Collapse time, s')
plt.show()
