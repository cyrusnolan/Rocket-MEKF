import numpy as np
from numpy.linalg import inv, norm
from scipy.linalg import block_diag
from consts import *
from utils import A, skew, Xi, rk4, jacob, qdot, wdot, rdot, vdot

##### Sim Data
sim_data = np.genfromtxt('sim_data.csv', delimiter=',')
sim_data = sim_data.T
fc_data = np.genfromtxt('fc_data.csv', delimiter=',')
fc_data = fc_data.T
time_arr = sim_data[0,:]
nk = len(time_arr)
# state vector
r_sim = sim_data[1:4,:]
q_sim = sim_data[4:8,:]
v_sim = sim_data[8:11,:]
w_sim = sim_data[11:14,:]
m_sim = sim_data[14:15,:]
X_sim = np.vstack([r_sim, q_sim, v_sim, w_sim, m_sim])
# control inputs
thrust = fc_data[7:8,:]
forces = sim_data[15:18,:]
moments = sim_data[18:21,:]

##### Create Measurements
# vector measurement 1
r_v1 = np.zeros([3,1])
r_v1[:,0] = (0, 0, 1)
sig_v1 = 0.01
R_v1 = sig_v1**2 * np.eye(3)
Z_v1 = np.zeros([3, nk])
Z_v1_clean = np.zeros_like(Z_v1)

# vector measurement 2
r_v2 = np.zeros([3,1])
r_v2[:,0] = (0, 1, 0)
sig_v2 = 0.005
R_v2 = sig_v2**2 * np.eye(3)
Z_v2 = np.zeros([3, nk])
Z_v2_clean = np.zeros_like(Z_v2)

for i in range(nk):
    q_true = q_sim[:, [i]]
    Z_v1_clean[:,[i]] = A(q_true)@r_v1
    Z_v1_unnorm = A(q_true)@r_v1 + np.random.normal(0, sig_v1, (3,1))
    Z_v1[:,[i]] = Z_v1_unnorm / norm(Z_v1_unnorm)
    Z_v2_clean[:,[i]] = A(q_true)@r_v2
    Z_v2_unnorm = A(q_true)@r_v2 + np.random.normal(0, sig_v2, (3,1))
    Z_v2[:,[i]] = Z_v2_unnorm / norm(Z_v2_unnorm)

# rate gyro
sig_rg = 0.003 # rad/s
R_rg = sig_rg**2 * np.eye(3)
Z_rg = w_sim + np.random.normal(0, sig_rg, w_sim.shape)

# gps
sig_gps = 10 # m
R_gps = sig_gps**2 * np.eye(3)
Z_gps = r_sim + np.random.normal(0, sig_gps, r_sim.shape)

# accel
sig_acc = 0.006 # m/s^2
R_acc = sig_acc**2 * np.eye(3)
Z_acc = forces + np.random.normal(0, sig_acc, forces.shape)

##### Initial Conditions
r0 = np.zeros([3,1])
q0 = np.zeros([4,1])
v0 = np.zeros([3,1])
w0 = np.zeros([3,1])
r0[:,0] = (-2500, 0, 2000)
q0[:,0] = (1, 0, 0, 0)
v0[:,0] = (50, 70, -75)
w0[:,0] = (0, 0, 0)
m0 = m_fuel + m_dry
P0 = np.diag([0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 100, 100, 100, 20, 20, 20])

##### run MEKF
# initialize states
nx = 12
# error rotation
dphi = np.zeros([3, nk])
# angular velocity
wu = np.zeros([3, nk])
wp = np.zeros_like(wu)
wp[:,[0]] = w0
# quaternion
qu = np.zeros([4, nk])
qp = np.zeros_like(qu)
qp[:,[0]] = q0
qnorm = np.zeros((nk,))
# position
ru = np.zeros([3,nk])
rp = np.zeros_like(ru)
rp[:,[0]] = r0
# velocity
vu = np.zeros([3,nk])
vp = np.zeros_like(vu)
vp[:,[0]] = v0
# covariance
Pu = np.zeros([nx, nx, nk])
Pp = np.zeros_like(Pu)
Pp[:,:,0] = P0

# initialize process noise
nw = 8
sig_w = 3000
Q = block_diag(sig_w**2 * np.eye(5), R_acc)
G = np.zeros([nx,nw])
G[3:5,0:2] = np.diag([1, 1])
G[9:12,5:8] = np.eye(3)

# initialize measurements
nz = 12
Z = np.vstack([Z_v1, Z_v2, Z_rg, Z_gps])
R = block_diag(R_v1, R_v2, R_rg, R_gps)
inn = np.zeros([nz,nk])
H_rg = np.zeros([3,nx])
H_rg[:,3:6] = np.eye(3)
H_gps = np.zeros([3,nx])
H_gps[:,6:9] = np.eye(3)

for k in range(nk):
    # gain
    b_v1 = A(qp[:,[k]])@r_v1
    H_v1 = np.zeros([3,nx])
    H_v1[:,:3] = skew(b_v1)
    b_v2 = A(qp[:,[k]])@r_v2
    H_v2 = np.zeros([3,nx])
    H_v2[:,:3] = skew(b_v2)
    H = np.vstack([H_v1, H_v2, H_rg, H_gps])
    K = Pp[:,:,k] @ H.T @ inv(H@Pp[:,:,k]@H.T + R)
    Kdphi = K[0:3,:]
    Kw = K[3:6,:]
    Kr = K[6:9,:]
    Kv = K[9:12,:]

    # update
    h = np.vstack([b_v1, b_v2, wp[:,[k]], rp[:,[k]]])
    inn[:,[k]] = Z[:,[k]] - h
    dphi[:,[k]] = Kdphi @ inn[:,[k]]
    wu[:,[k]] = wp[:,[k]] + Kw@inn[:,[k]]
    ru[:,[k]] = rp[:,[k]] + Kr@inn[:,[k]]
    vu[:,[k]] = vp[:,[k]] + Kv@inn[:,[k]]
    Pu[:,:,k] = (np.eye(nx) - K@H) @ Pp[:,:,k]

    # reset
    qu_unnorm = qp[:,[k]] + 1/2*Xi(qp[:,[k]])@dphi[:,[k]]
    qnorm[k] = norm(qu_unnorm)
    qu[:,[k]] = qu_unnorm / qnorm[k]

    # propagate
    if k < nk-1:
        thrust_idx = int(np.floor(k/10))
        qp[:,[k+1]] = rk4(qdot, qu[:,[k]], dt, wu[:,[k]])
        wp[:,[k+1]] = rk4(wdot, wu[:,[k]], dt, moments)
        rp[:,[k+1]] = rk4(rdot, ru[:,[k]], dt, vu[:,[k]])
        vp[:,[k+1]] = rk4(vdot, vu[:,[k]], dt, Z_acc[:,[k]], m_sim[0,k])
        Fc = jacob(np.vstack([dphi[:,[k]], wp[:,[k]]]))
        Fd = rk4(lambda X, Fc: Fc @ X, np.eye(nx), dt, Fc)
        G[9:12,2:5] = np.eye(3) / m_sim[0,k]
        Qk = G@Q@G.T * dt
        Pp[:,:,k+1] = Fd@Pu[:,:,k]@Fd.T + Qk
