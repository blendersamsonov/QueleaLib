import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
# motion.py - интерфейс к нашей библиотеке на Rust
import motion

lw = 0.3 # linewidth
lwl = 1.2*lw # linewidth for lines in plots

fig = plt.figure(figsize = (3.3, 3.3))

font = {'family' : 'Liberation Serif',
        'weight' : 'normal',
	'size'   : 9}
matplotlib.rc('font', **font)
matplotlib.rc('lines', linewidth=lw)
matplotlib.rc('axes', linewidth=lw)

proj = '3d'
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111, projection=proj)
ax.view_init(elev = 25, azim = 45-180)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$z$')

ax.w_xaxis.set_pane_color((1, 1, 1, 1))
ax.w_yaxis.set_pane_color((1, 1, 1, 1))
ax.w_zaxis.set_pane_color((1, 1, 1, 1))

def save_plot(dots, labels, limits, elev, azim, filename, type):
        lw = 0.3
        lwl = 1.2*lw

        fig = plt.figure(figsize = (3.3, 3.3))

        font = {'family' : 'Liberation Serif',
                'weight' : 'normal',
                'size'   : 9}
        matplotlib.rc('font', **font)
        matplotlib.rc('lines', linewidth=lw)
        matplotlib.rc('axes', linewidth=lw)

        ax = fig.add_subplot(111, projection=proj)
        ax.view_init(elev = elev, azim = azim-180)

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])

        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.set_zlim(limits[4], limits[5])

        ax.w_xaxis.set_pane_color((1, 1, 1, 1))
        ax.w_yaxis.set_pane_color((1, 1, 1, 1))
        ax.w_zaxis.set_pane_color((1, 1, 1, 1))

        if type == 0:
                plt.plot(dots[:,0], dots[:,1], dots[:,2], '-b', linewidth = lwl)
        if type == 1:
                plt.plot(dots[:,3], dots[:,4], dots[:,5], '-b', linewidth = lwl)
        if type == 2:
                plt.plot(dots[:,0], dots[:,1], dots[:,2], '-b', linewidth = lwl)
                plt.plot(dots[:,6], dots[:,7], dots[:,8], '--r', linewidth = lwl)
        if type == 3:
                plt.plot(dots[:,9], dots[:,10], 0, '-b', linewidth = lwl)
        plt.savefig(filename+'.pdf')
        plt.close()

def piazza_solution(type, omega, phi0, E0, v0, dt, N, coords):
        gamma0 = 1 / math.sqrt(1 - v0 * v0)
        xi = E0 / omega
        eta = omega * gamma0 * (1 + v0)
        if type == 'linear':
                def H(phi):
                        cos = math.cos(phi)
                        sin = math.sin(phi)
                        cos0 = math.cos(phi0)
                        sin0 = math.sin(phi0)
                        return 1 + xi*xi*eta*(phi - phi0 - sin * cos + sin0 * cos0)/(3*137)
                
                def I1(phi):
                        h = H(phi)
                        cos = math.cos(phi)
                        sin = math.sin(phi)
                        cos0 = math.cos(phi0)
                        sin0 = math.sin(phi0)
                        return cos * h - cos0 - 2 * eta * (1 + xi*xi*(sin*sin + sin*sin0 + sin0*sin0)/3)*(sin - sin0)/(3*137)
                
                def V(phi):
                        h = H(phi)
                        i1 = I1(phi)
                        tmp = omega * (h * h - 1 + xi*xi*i1*i1)/(2*eta)
                        return [gamma0 + tmp,
                                0,
                                -v0*gamma0 + tmp,
                                -xi*i1]
        
                path = np.empty([N+1,3])
                path[0] = [coords[1],coords[2],coords[3]]
                for i in range(N):
                        phi = phi0 + omega * (coords[0] - coords[2]) #omega(t-y)
                        h = H(phi)
                        v = V(phi)
                        for j in range(4):
                                coords[j] = coords[j] + h * v[j] * dt / (gamma0*(1+v0))
                                if j != 0:
                                        path[N, j-1] = coords[j]
                return path
        

pi = math.pi
#параметры для drive_particle: 
#нач. координаты, нач. скорость, 
#направление, амплитуда, частота, волновое число, напр-е волнового вектора, эллипс поляризации [по полю, перп. полю], сдвиг фазы E и B (задаётся на правой границе импульса)
#огибающая (0 - константа, 1 - косинус, 2 - нет границ), размеры импульса [x, y, z], где x - по напр-ю k (!!!пока по глобальным осям!!!), центр импульса
#dt, число шагов, схема рад. трения
#c = 1. Для э-м волны omega = k

test = 2
###ТЕСТЫ:
#
#1. Постоянное магнитное поле (сильный релятивизм (gamma = 707), без трения)
if test == 1:
        limits = [-10, 10,
                  -10, 10,
                  -1, 1]
        p_limits = [-10, 10,
                  -10, 10,
                  -1, 1]

        v = motion.drive_particle([-6, 6, 0], [0.999999, 0.0, 0.0],
                                [0, 0, 0], 0, 0, 0, [1, 0, 0], [1, 0], 0,
                                [0, 0, 1], 50, 0, 0, [1, 0, 0], [1, 0], 0,
                                1, [10000, 10000, 10000], [0, 0, 0],
                                0.0001, 10000000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)


#2. Постоянное магнитное поле (сильный релятивизм, трение)
if test == 2:
        limits = [-7, -3,
                  3, 7,
                  -1, 5]
        p_limits = [-10, 10,
                  -10, 10,
                  -1, 1]

        #steps = 10
        steps = 1000000
        dt = 0.0001
        vel = [0.9999999, 0.0, 0.0]
        gamma0 = 1/math.sqrt(1-vel[0]*vel[0]-vel[1]*vel[1]-vel[2]*vel[2])

        gamma_limits=[1,gamma0,
                      0,dt*(steps+1),
                      -1,1]

        B = 500
        v = motion.drive_particle([-6, 6, 0], vel,
                                [0, 0, 0], 0, 0, 0, [1, 0, 0], [1, 0], 0,
                                [0, 0, 1], B, 0, 0, [1, 0, 0], [1, 0], 0,
                                1, [10000, 10000, 10000], [0, 0, 0],
                                dt, steps, 1)
        gamma_th = np.empty([steps+1,3])
        gamma_th[0]=[gamma0,0,0]
        for i in range(0,steps):
                a = math.tanh(B*B*dt*dt*dt*i/6) # поле B*dt/2 (как в схеме Вэя)
                gamma_th[i+1] = [(gamma0+a)/(1+gamma0*a),dt*i,0]
        #print(gamma_th)
        gamma_dif = np.empty([steps,3])
        for i in range(0,steps):
                gamma_dif[i] = [gamma_th[i,0] - v[i,9],dt*i,0]
        save_plot(v, ['x','y','z'], limits, 90, 0, 'particle_track', 0)
        save_plot(gamma_th,['gamma','t',''], gamma_limits, 90, 0,'gamma_th',0)
        save_plot(gamma_dif,['gamma','t',''], [-10,10,0,dt*steps,-1,1], 90, 0,'gamma_dif',0)
        save_plot(v,['gamma','t',''], gamma_limits, 90, 0,'gamma', 3)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#3. Радиоимпульс электрического поля (огибающая - cos^2(x^4)). !!!Странный дрейф по z. С огибающей - прямоугольником такого нет.
if test == 3:
        limits = [-10, 0,
                  -10, 10,
                  -1, 1]
        p_limits = [-10, 0,
                  -10, 10,
                  -1, 1]
        
        omega = 1
        v = motion.drive_particle([0, 0, 0], [-0.1, 0, 0],
                                  [0, 0, 1], 1, omega, omega, [1, 0, 0], [1, 0], 0,
                                  [0, 1, 0], 0, 0, 0, [1, 0, 0], [1, 0], 0,
                                  1, [100.0, 1000.0, 1000.0], [-55.0, 0.0, 0.0],
                                  0.0001, 1000000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#4. Линейно-поляризованная электромагнитная волна (огибающая - cos^2(x^4))
if test == 4:
        limits = [-10, 0,
                  -10, 10,
                  -0.7, 0.7]
        p_limits = [-5, 0,
                  0, 5,
                  0, 5]

        
        omega = 1.5
        v = motion.drive_particle([0, 0, 0], [0.0, 0, 0],
                                  [0, 0, 1], 20, omega, omega, [1, 0, 0], [1, 0], 0,
                                  [0, 1, 0], 20, omega, omega, [1, 0, 0], [1, 0], 0,
                                  1, [120.0, 1000.0, 1000.0], [-60.0, 0.0, 0.0],
                                  0.0001, 1500000, 1)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        save_plot(v, ['px','py','pz'], p_limits, 0, 90, 'particle_impulse', 1)

#5. Линейно-поляризованная электромагнитная волна (огибающая - прямоугольник, частица "убегает" от волны)
if test == 5:
        limits = [-4, 0,
                  -10, 10,
                  -0.7, 0.7]
        p_limits = [-4, 0,
                  -10, 10,
                  -0.7, 0.7]
        
        omega = 1.5
        v = motion.drive_particle([0, 0, 0], [0.06, 0, 0],
                                  [0, 0, 1], 1, omega, omega, [1, 0, 0], [1, 0], 0,
                                  [0, 1, 0], 1, omega, omega, [1, 0, 0], [1, 0], 0,
                                  0, [120.0, 1000.0, 1000.0], [-60.0, 0.0, 0.0],
                                  0.0001, 1000000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#6. Циркулярно-поляризованная электромагнитная волна (огибающая - косинус). !!!Непостоянный радиус или центр спирали
if test == 6:
        limits = [-8, 0,
                  -1, 1,
                  -1, 1]
        p_limits = [-8, 0,
                  -1, 1,
                  -1, 1]
        
        omega = 1.5
        v = motion.drive_particle([0, 0, 0], [0.00, 0, 0],
                                  [0, 0, 1], 1, omega, omega, [1, 0, 0], [1, 1], 0,
                                  [0, 1, 0], 1, omega, omega, [1, 0, 0], [1, 1], 0,
                                  1, [120.0, 1000.0, 1000.0], [-60.0, 0.0, 0.0],
                                  0.0001, 1100000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#7. Циркулярно-поляризованная электромагнитная волна (огибающая - косинус). Сильный релятивизм. 
if test == 7:
        limits = [-15, 0,
                  -0.3, 0.3,
                  -0.3, 0.3]
        p_limits = [-15, 0,
                  -0.3, 0.3,
                  -0.3, 0.3]
        
        omega = 3
        v = motion.drive_particle([0, 0, 0], [0.00, 0, 0],
                                  [0, 0, 1], 10000, omega, omega, [1, 0, 0], [1, 1], 0,
                                  [0, 1, 0], 10000, omega, omega, [1, 0, 0], [1, 1], 0,
                                  1, [120.0, 1000.0, 1000.0], [-60.0, 0.0, 0.0],
                                  0.0001, 1100000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#8. Циркулярно-поляризованная электромагнитная волна (огибающая - косинус). Сильный релятивизм. Трение
if test == 8:
        limits = [-15, 0,
                  -0.3, 0.3,
                  -0.3, 0.3]
        p_limits = [-15, 0,
                  -0.3, 0.3,
                  -0.3, 0.3]
        
        omega = 3
        v = motion.drive_particle([0, 0, 0], [0.00, 0, 0],
                                  [0, 0, 1], 10000, omega, omega, [1, 0, 0], [1, 1], 0,
                                  [0, 1, 0], 10000, omega, omega, [1, 0, 0], [1, 1], 0,
                                  1, [120.0, 1000.0, 1000.0], [-60.0, 0.0, 0.0],
                                  0.0001, 1100000, 1)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 0)
        #save_plot(v, ['px','py','pz'], p_limits, 25, 45, 'particle_impulse', 1)

#9. Линейно-поляризованная электромагнитная волна (огибающая - степ функция). Сильный релятивизм. Трение
if test == 9:
        limits = [-20, 0,
                  -10, 10,
                  0, 1]
        p_limits = [-5, 0,
                  0, 5,
                  0, 5]

        
        omega = 1
        v = motion.drive_particle([0, 0, 0], [0.0, 0, 0],
                                  [0, 0, 1], 1000, omega, omega, [1, 0, 0], [1, 0], 0,
                                  [0, 1, 0], 1000, omega, omega, [1, 0, 0], [1, 0], 0,
                                  2, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                                  0.00001, 5000000, 0)
        save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track', 2)
        save_plot(v, ['px','py','pz'], p_limits, 0, 90, 'particle_impulse', 1)
        #save_plot(v, ['x','y','z'], limits, 25, 45, 'particle_track_stationary', 2)


        #v = piazza_solution(type = 'linear', omega = omega, phi0 = 0, E0 = 10000, v0 = 0, dt = 0.01, N = 100, coords = [0,0,0,0])
        #print(v)
        #save_plot(v, ['x','y','z'], limits, 0, 0, 'di_piazza_solution', 0)