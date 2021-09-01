import sys
sys.path.append("../figs/")
import motion

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

T = 200
dt = 0.005

_lambda1 = 0.5
a01 = 5000
_1_k1 = _lambda1 / (2 * np.pi)
coef1 = np.power(_1_k1 / 2.82e-9, 1/3)

_lambda2 = 1000
a02 = 1500
_1_k2 = _lambda2 / (2 * np.pi)
coef2 = np.power(_1_k2 / 2.82e-9, 1/3)

print(coef1, a01)
print(coef2, a02)

dp = motion.drive_particle
traj1 = dp([0,0,0,0,0,0], 0, dt, int(T / dt), 3, [a01, 1, 0, 0, 0, 0, 1], _lambda1)
traj1 = pd.DataFrame(traj1, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
traj1.g = np.sqrt(1 + traj1.ux**2 + traj1.uy**2 + traj1.uz**2)
traj1.t = traj1.index * dt
traj1.ez = -np.cos(traj1.t)
traj1.ey = np.sin(traj1.t)
traj1.p = np.sqrt(traj1.uy**2+traj1.uz**2)
traj1.pz = traj1.uz / traj1.p
traj1.py = traj1.uy / traj1.p
traj1.theta = 2*np.arccos(traj1.ez * traj1.pz + traj1.ey * traj1.py)/np.pi

traj2 = dp([0,0,0,0,0,0], 0, dt, int(T / dt), 3, [a02, 1, 0, 0, 0, 0, 1], _lambda2)
traj2 = pd.DataFrame(traj2, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
traj2.g = np.sqrt(1 + traj2.ux**2 + traj2.uy**2 + traj2.uz**2)
traj2.t = traj2.index * dt
traj2.ez = -np.cos(traj2.t)
traj2.ey = np.sin(traj2.t)
traj2.p = np.sqrt(traj2.uy**2+traj2.uz**2)
traj2.pz = traj2.uz / traj2.p
traj2.py = traj2.uy / traj2.p
traj2.theta = 2*np.arccos(traj2.ez * traj2.pz + traj2.ey * traj2.py)/np.pi


lw = 0.5 # linewidth
lwl = 1.2*lw # linewidth for lines in plots

fig = plt.figure(figsize = (5, 10))

font = {'family' : 'Liberation Serif',
        'weight' : 'normal',
	'size'   : 9}
matplotlib.rc('font', **font)
matplotlib.rc('lines', linewidth=lw)
matplotlib.rc('axes', linewidth=lw)

ax = fig.add_subplot(211)

#ax.set_xlim(-2, 2)
#ax.set_ylim(-2, 2)
#ax.set_zlim(-1.5, 1.5)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

#ax.w_xaxis.set_pane_color((1, 1, 1, 1))
#ax.w_yaxis.set_pane_color((1, 1, 1, 1))

plt.plot(traj1.z, traj1.y, '-r', linewidth = lwl)
plt.plot(traj2.z, traj2.y, '-g', linewidth = lwl)

#ax2 = fig.add_subplot(312)
#plt.plot(traj1.t, traj1.g, '-r', linewidth = lwl)
#plt.xlabel(r'$t$')
#plt.ylabel(r'$g$')

ax3 = fig.add_subplot(212)
plt.plot(traj1.t, traj1.theta, '-r', linewidth = lwl)
plt.plot(traj2.t, traj2.theta, '-g', linewidth = lwl)
ax3.set_ylim(0, 1)
plt.xlabel(r'$t$')
plt.ylabel(r'$2phi/pi$')

plt.savefig("rot_E.pdf")
