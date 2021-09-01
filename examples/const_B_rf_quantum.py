import sys
sys.path.append("../figs/")
import motion

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

# Аналитическое решение для энергии в приближении синхротронного излучения:
# g = (1 + g0 * cth(A * t)) / (g0 + cth(A * t)), где
# A = 4 * Pi * re * a0^2 / (3 * lambda), где 
# re - классический радиус электрона: 2.82*10^-13 cm
# lambda - длина волны, на которую нормируем (в коде 1 um, если не задавать)
# время t нормированно на частоту omega = 2 * Pi * c / lambda
# Для совпадения с теорией и отсутствия скручивания в режиме без трения необходим достаточно малый шаг dt

T = 200
dt = 0.005
_lambda = 0.1
a0 = 50
v = 0.99999
g0 = 1/np.sqrt(1-v*v)
p = v * g0
A = 4 * np.pi * 2.82e-13 * a0 * a0 / (3 * _lambda * 1e-4)

dp = motion.drive_particle
traj1 = dp([0,0,0,p,0,0], dt, int(T / dt), 2, [a0], _lambda, 0)
traj1 = pd.DataFrame(traj1, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
traj1.g = np.sqrt(1 + traj1.ux**2 + traj1.uy**2 + traj1.uz**2)
traj1.t = traj1.index * dt

traj2 = dp([0,0,0,p,0,0], dt, int(T / dt), 2, [a0], _lambda, 1)
traj2 = pd.DataFrame(traj2, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
traj2.g = np.sqrt(1 + traj2.ux**2 + traj2.uy**2 + traj2.uz**2)
traj2.t = traj2.index * dt

traj_av = dp([0,0,0,p,0,0], dt, int(T / dt), 2, [a0], _lambda, 1)
traj_av = pd.DataFrame(traj_av, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
traj_av.g = np.sqrt(1 + traj_av.ux**2 + traj_av.uy**2 + traj_av.uz**2)
steps = 200
for i in range(1,steps-1):
        traj3 = dp([0,0,0,p,0,0], dt, int(T / dt), 2, [a0], _lambda, 1)
        traj3 = pd.DataFrame(traj3, columns = ['x', 'y', 'z', 'ux', 'uy', 'uz'])
        traj3.g = np.sqrt(1 + traj3.ux**2 + traj3.uy**2 + traj3.uz**2)
        traj_av.g = traj_av.g + traj3.g
        traj_av.x = traj_av.x + traj3.x
        traj_av.y = traj_av.y + traj3.y
traj_av.g = traj_av.g / steps
traj_av.x = traj_av.x / steps
traj_av.y = traj_av.y / steps

print("g_amp = ", max(traj2.g))

lw = 0.5 # linewidth
lwl = 1.2*lw # linewidth for lines in plots

fig = plt.figure(figsize = (5, 5))

font = {'family' : 'Liberation Serif',
        'weight' : 'normal',
	'size'   : 9}
matplotlib.rc('font', **font)
matplotlib.rc('lines', linewidth=lw)
matplotlib.rc('axes', linewidth=lw)

#ax = fig.add_subplot(111)

#ax.set_xlim(-2, 2)
#ax.set_ylim(-2, 2)
#ax.set_zlim(-1.5, 1.5)

#ax.set_xlabel(r'$x$')
#ax.set_ylabel(r'$y$')

#ax.w_xaxis.set_pane_color((1, 1, 1, 1))
#ax.w_yaxis.set_pane_color((1, 1, 1, 1))

#plt.plot(traj1.x, traj1.y, '--r', linewidth = 0.1*lwl)
#plt.plot(traj2.x, traj2.y, '-g', linewidth = 0.5*lwl)
#plt.plot(traj_av.x, traj_av.y, '-b', linewidth = 0.1*lwl)

ax2 = fig.add_subplot(111)
#ax2.set_ylim(1, g0+10)

plt.plot(traj1.t, traj1.g, '-r', linewidth = lwl)
plt.plot(traj2.t, traj2.g, '-g', linewidth = lwl)
plt.plot(traj2.t, traj_av.g, '-b', linewidth = lwl)

plt.xlabel(r'$t$')
plt.ylabel(r'$g$')

plt.savefig("const_B_q.pdf")
