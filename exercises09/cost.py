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
plt.show()
