import numpy as np
from consts import *
from typing import Callable

def eom(X: np.ndarray, F: np.ndarray, M: np.ndarray, T: float):
    # r = X[0:3]
    q = X[3:7]
    v = X[7:10]
    w = X[10:13]
    m = X[13]

    rdot = v
    qdot = 0.5 * crassidis(q,w)
    vdot = F / m
    w1dot = 1 / I1 * (M[0,0] + w[1,0]*w[2,0]*(I2 - I3))
    w2dot = 1 / I2 * (M[1,0] + w[2,0]*w[0,0]*(I3 - I1))
    w3dot = 1 / I3 * (M[2,0] + w[0,0]*w[1,0]*(I1 - I2))
    wdot = np.vstack([w1dot,w2dot,w3dot])
    mdot = -T / (g0 * Isp)
    return np.vstack([rdot,qdot,vdot,wdot,mdot])

def hamilton(q: np.ndarray, r: np.ndarray):
    p0 = r[0,0]*q[0,0] - r[1,0]*q[1,0] - r[2,0]*q[2,0] - r[3,0]*q[3,0]
    p1 = r[0,0]*q[1,0] + r[1,0]*q[0,0] - r[2,0]*q[3,0] + r[3,0]*q[2,0]
    p2 = r[0,0]*q[2,0] + r[1,0]*q[3,0] + r[2,0]*q[0,0] - r[3,0]*q[1,0]
    p3 = r[0,0]*q[3,0] - r[1,0]*q[2,0] + r[2,0]*q[1,0] + r[3,0]*q[0,0]
    return np.vstack([p0,p1,p2,p3])

def crassidis(q: np.ndarray, r:np.ndarray):
    p0 = -r[0,0]*q[1,0] - r[1,0]*q[2,0] - r[2,0]*q[3,0]
    p1 =  r[0,0]*q[0,0] - r[1,0]*q[3,0] + r[2,0]*q[2,0]
    p2 =  r[0,0]*q[3,0] + r[1,0]*q[0,0] - r[2,0]*q[1,0]
    p3 = -r[0,0]*q[2,0] + r[1,0]*q[1,0] + r[2,0]*q[0,0]
    return np.vstack([p0,p1,p2,p3])

def rk4(X: np.ndarray, F: np.ndarray, M: np.ndarray, T: float, h: float,
        f: Callable[[np.ndarray, np.ndarray, np.ndarray, float], np.ndarray]):

    k1 = f(X, F, M, T)              # f(t, X)
    k2 = f(X + h*k1/2, F, M, T)     # f(t + h/2, X + h*k1/2)
    k3 = f(X + h*k2/2, F, M, T)     # f(t + h/2, X + h*k2/2)
    k4 = f(X + h*k3, F, M, T)       # f(t + h, X + h*k3)
    return X + h/6*(k1 + 2*k2 + 2*k3 +k4)
