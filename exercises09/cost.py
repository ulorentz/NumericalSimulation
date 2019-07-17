#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
ave=np.loadtxt("av_circle.dat")
best=np.loadtxt("best_circle.dat")

fig,ax = plt.subplots(1, 2, figsize=(16, 6))

ax[0].loglog(ave, label="Averaged")
ax[0].loglog(best, label="Best")
ax[0].set_title("circle")
ax[0].legend()
ave=np.loadtxt("av_squared.dat")
best=np.loadtxt("best_squared.dat")
ax[1].loglog(ave, label="Averaged")
ax[1].loglog(best, label="Best")
ax[1].set_title("squared")
ax[1].legend();


fig,ax = plt.subplots(1, 1, figsize=(6, 6))

x, y=np.loadtxt("outputs/squared_cities_location.dat", delimiter="\t", unpack=True)
final=np.loadtxt("outputs/squared_final_best_path.dat",dtype="uint32")  
x_final=x[final]
x_final=np.append(x_final, x[final[0]])
y_final=y[final]
y_final=np.append(y_final, y[final[0]])
ax.plot(x_final, y_final, marker=".");  
ax.set_title("Final best L2 path");


plt.show()
