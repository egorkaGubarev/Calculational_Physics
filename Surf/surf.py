import numpy as np
import matplotlib.pyplot as plt

draw_connect = True

circ = np.loadtxt("../../../VSProjects/Calculational_physics/Surf_no_refer/circle.txt")

dots_on_circ = int(circ[0, 0])
soft_nodes = int(circ[0, 1])
soft_string_number = dots_on_circ + 1
connect_string_number = dots_on_circ + soft_nodes + 1

circ_x = circ[1:soft_string_number, 0]
circ_y = circ[1:soft_string_number, 1]
circ_z = circ[1:soft_string_number, 2]

# Close circuit
circ_x = np.append(circ_x, circ[1, 0])
circ_y = np.append(circ_y, circ[1, 1])
circ_z = np.append(circ_z, circ[1, 2])

nodes_x = circ[soft_string_number: connect_string_number, 0]
nodes_y = circ[soft_string_number: connect_string_number, 1]
nodes_z = circ[soft_string_number: connect_string_number, 2]

connect_x = circ[connect_string_number:, 0]
connect_y = circ[connect_string_number:, 1]
connect_z = circ[connect_string_number:, 2]

connect_strings_numb = len(connect_x)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_title('Circuit')
ax.plot(circ_x, circ_y, circ_z)
ax.plot(nodes_x, nodes_y, nodes_z, '.')

if draw_connect:
    for i in range(0, connect_strings_numb - 1, 2):
        connect_index = i + 1

        node_x = connect_x[i]
        node_y = connect_y[i]
        node_z = connect_z[i]

        connected_x = connect_x[connect_index]
        connected_y = connect_y[connect_index]
        connected_z = connect_z[connect_index]

        x = (node_x, connected_x)
        y = (node_y, connected_y)
        z = (node_z, connected_z)

        ax.plot(x, y, z, color='red')

plt.show()
