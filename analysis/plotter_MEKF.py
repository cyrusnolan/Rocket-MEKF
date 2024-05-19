import numpy as np
import matplotlib.pyplot as plt
from MEKF import time_arr
from MEKF import Z
from MEKF import dphi, Pu, inn, q_sim, qu, qnorm, wu, Z_rg, w_sim, ru, Z_gps, r_sim, vu, v_sim

# # dphi
# for i in range(3):
#     sig = 2*np.sqrt(Pu[i,i,:])
#     plt.figure()
#     plt.plot(time_arr, dphi[i,:])
#     plt.fill_between(time_arr, 2*sig, -2*sig, color='gray', alpha=0.2, label='Error bound')
# plt.show()

# # autocorrelation of the innovations
# plt.figure()
# plt.acorr(inn[0,:], maxlags=30, normed=True)
# plt.title('Autocorrelation')

# # quaternion estimate
# for i in range(4):
#     plt.figure()
#     plt.plot(time_arr, q_sim[i, :], label='q_sim')
#     plt.plot(time_arr, qu[i, :], label='qu')
#     plt.xlabel('Time')
#     plt.ylabel('First Element')
#     plt.legend()
# plt.show()

# # q norm
# plt.figure()
# plt.plot(time_arr, qnorm)
# ax = plt.gca()
# ax.get_yaxis().get_major_formatter().set_useOffset(False)
# plt.show()

# # angular velocity
# for i in range(3):
#     sig = np.sqrt(Pu[3+i, 3+i, :])
#     plt.figure()
#     plt.plot(time_arr, wu[i,:], label='estimate')
#     plt.plot(time_arr, Z_rg[i,:], label='measurement')
#     plt.plot(time_arr, w_sim[i,:], label='true value')
#     upper_bound = wu[i,:] + 2*sig
#     lower_bound = wu[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.legend()
# plt.show()

# angular veloctiy error
for i in range(3):
    sig = np.sqrt(Pu[3+i, 3+i, :])
    err = wu[i,:] - w_sim[i,:]
    plt.figure()
    plt.plot(time_arr, err, label='estimate error')
    upper_bound = 2*sig
    lower_bound = -2*sig
    plt.fill_between(time_arr, lower_bound, upper_bound,
                     color='gray', alpha=0.2, label='$\pm 2\sigma$')
    plt.legend()
plt.show()

# # position
# for i in range(3):
#     sig = np.sqrt(Pu[6+i,6+i,:])
#     plt.figure()
#     plt.plot(time_arr, ru[i,:], label='estimate')
#     plt.plot(time_arr, Z_gps[i,:], label='measurement')
#     plt.plot(time_arr, r_sim[i,:], label='true value')
#     upper_bound = ru[i,:] + 2*sig
#     lower_bound = ru[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.legend()
# plt.show()

# # position error
# for i in range(3):
#     sig = np.sqrt(Pu[6+i,6+i,:])
#     err = ru[i,:] - r_sim[i,:]
#     plt.figure()
#     plt.plot(time_arr, err, label='estimate error')
#     upper_bound = 2*sig
#     lower_bound = -2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.legend()
# plt.show()

# # velocity
# for i in range(3):
#     sig = np.sqrt(Pu[9+i,9+i,:])
#     plt.figure()
#     plt.plot(time_arr, vu[i,:], label='estimate')
#     plt.plot(time_arr, v_sim[i,:], label='true value')
#     upper_bound = vu[i,:] + 2*sig
#     lower_bound = vu[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound, 
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.legend()
# plt.show()

# # velocity error
# for i in range(3):
#     sig = np.sqrt(Pu[9+i,9+i,:])
#     err = vu[i,:] - v_sim[i,:]
#     plt.figure()
#     plt.plot(time_arr, err, label='estimate error')
#     upper_bound = 2*sig
#     lower_bound = -2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.legend()
# plt.show()
