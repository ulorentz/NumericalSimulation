#!/usr/bin/python
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
x100,y100,z100 = np.loadtxt("100.txt", delimiter="\t", unpack=True)
x210,y210,z210 = np.loadtxt("210.txt", delimiter="\t", unpack=True)


fig = plt.figure(figsize=(10,8))
ax = Axes3D(fig)
ax.scatter(x100, y100, z100, c=z100, marker='.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.view_init(10, 30)

fig = plt.figure(figsize=(10,8))
ax = Axes3D(fig)
ax.scatter(x210, y210, z210, c=z210, marker='.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.view_init(10, 30)

plt.show()
