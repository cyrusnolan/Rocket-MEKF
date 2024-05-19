import numpy as np
import matplotlib.pyplot as plt
from consts import *
from utils import rk4, eom


# import control inputs
sim_data = np.genfromtxt('sim_data.csv', delimiter=',')
fc_data = np.genfromtxt('fc_data.csv', delimiter=',')
t_arr = sim_data[:,0]
thrust = fc_data[:,7:8]
forces = sim_data[:,15:18]
moments = sim_data[:,18:21]

# Initial Conditions
r0 = np.zeros((3,1))
q0 = np.zeros((4,1))
v0 = np.zeros((3,1))
w0 = np.zeros((3,1))
r0[:,0] = (-2500, 0, 2000)
q0[:,0] = (1, 0, 0, 0)
v0[:,0] = (50, 70, -75)
w0[:,0] = (0, 0, 0)
m0 = m_fuel + m_dry
X0 = np.vstack([r0, q0, v0, w0, m0])

# run sim
nx = len(X0)
nt = len(t_arr)
X = np.zeros((nx, nt))
X[:,0:1] = X0
for n, t in enumerate(t_arr):
    k = int(np.floor(n / freq))
    T = thrust[k].item()
    F = forces[n,:].reshape(3,1)
    M = moments[n,:].reshape(3,1)
    X[:,n+1:n+2] = rk4(X[:,n:n+1], F, M, T, dt, eom)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(X[0, :], X[1, :], X[2, :])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.grid(True)
ax.set_box_aspect([1, 1, 1])
plt.show()
