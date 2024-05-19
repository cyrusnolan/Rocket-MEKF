import numpy as np
from consts import *
from typing import Any, Callable

def eom(X: np.ndarray, Force: np.ndarray, M: np.ndarray, Thrust: float):
    # r = X[0:3]
    q = X[3:7]
    v = X[7:10]
    w = X[10:13]
    m = X[13]

    rdot = v
    qdot = 0.5 * Xi(q)@w
    vdot = Force / m
    w1dot = 1 / I1 * (M[0,0] + w[1,0]*w[2,0]*(I2 - I3))
    w2dot = 1 / I2 * (M[1,0] + w[2,0]*w[0,0]*(I3 - I1))
    w3dot = 1 / I3 * (M[2,0] + w[0,0]*w[1,0]*(I1 - I2))
    wdot = np.vstack([w1dot,w2dot,w3dot])
    mdot = -Thrust / (g0 * Isp)
    return np.vstack([rdot,qdot,vdot,wdot,mdot])

def qdot(q: np.ndarray, w: np.ndarray):
    return 0.5 * Xi(q)@w

def wdot(w: np.ndarray, M: np.ndarray):
    w1dot = 1 / I1 * (M[0,0] + w[1,0]*w[2,0]*(I2 - I3))
    w2dot = 1 / I2 * (M[1,0] + w[2,0]*w[0,0]*(I3 - I1))
    w3dot = 1 / I3 * (M[2,0] + w[0,0]*w[1,0]*(I1 - I2))
    return np.vstack([w1dot,w2dot,w3dot])

def rdot(r: np.ndarray, v: np.ndarray):
    return v

def vdot(v: np.ndarray, z_acc: np.ndarray, m: float):
    return z_acc/m

def Xi(q: np.ndarray):
    return np.vstack([(-q[1:].T), (q[0,0]*np.eye(3) + skew(q[1:]))])

def skew(vv: np.ndarray):
    v = vv.squeeze()
    return np.array([(0, -v[2], v[1]), (v[2], 0, -v[0]), (-v[1], v[0], 0)])

def A(q: np.ndarray):
    return (q[0,0]**2 - (q[1,0]**2 + q[2,0]**2 + q[3,0]**2))*np.eye(3) - 2*q[0,0]*skew(q[1:]) + 2*q[1:]@q[1:].T

def jacob(x):
    dphi = x[0:3]
    w = x[3:6]
    J = np.zeros((12,12))
    J[0:3,0:3] = -skew(w)
    J[0:3,3:6] = skew(dphi)
    J[3,3:6] = (0, (I2 - I3) / I1 * w[2,0], (I2 - I3) / I1 * w[1,0])
    J[4,3:6] = ((I3 - I1) / I2 * w[2,0], 0, (I3 - I1) / I2 * w[0,0])
    J[5,3:6] = ((I1 - I2) / I3 * w[1,0], (I1 - I2) / I3 * w[0,0], 0)
    J[6:9,9:12] = np.eye(3)
    return J

def rk4(f: Callable[[np.ndarray, Any], np.ndarray], X: np.ndarray, dt: float, *args):
    k1 = f(X, *args)                # f(t, X)
    k2 = f(X + dt*k1/2, *args)      # f(t + h/2, X + h*k1/2)
    k3 = f(X + dt*k2/2, *args)      # f(t + h/2, X + h*k2/2)
    k4 = f(X + dt*k3, *args)        # f(t + h, X + h*k3)
    return X + dt/6*(k1 + 2*k2 + 2*k3 +k4)
