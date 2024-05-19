import numpy as np
import matplotlib.pyplot as plt
from MEKF import time_arr
from MEKF import Z
from MEKF import dphi, Pu, inn, q_sim, qu, qnorm, wu, Z_rg, w_sim, ru, Z_gps, r_sim, vu, v_sim
from scipy.signal import welch
plt.rcParams['font.size'] = 20

# # dphi
# for i in range(3):
#     sig = 2*np.sqrt(Pu[i,i,:])
#     plt.figure()
#     plt.plot(time_arr, dphi[i,:])
#     plt.fill_between(time_arr, 2*sig, -2*sig, color='gray', alpha=0.2, label='Error bound')

# autocorrelation of the innovations
for i in range(3):
    lab = ['$\delta \phi_1$', '$\delta \phi_2$', '$\delta \phi_3$']
    plt.figure()
    plt.acorr(inn[i,:], maxlags=30, normed=True)
    plt.title(f'Autocorrelation: {lab[i]}')
    plt.xlabel("Time Steps, Frequecy = 0.1 Hz")
plt.show()

# # Power Spectrum Test
# frequencies, psd = welch(inn[0,:])
# plt.figure()
# plt.semilogy(frequencies, psd)
# plt.title('Power spectral density')
# plt.xlabel('Frequency')
# plt.ylabel('PSD')
# plt.show()

# # quaternion estimate
# for i in range(4):
#     lab = ['q1 (scalar part)', 'q2', 'q3', 'q4']
#     plt.figure()
#     plt.plot(time_arr, qu[i, :], label='estimate')
#     plt.plot(time_arr, q_sim[i, :], label='true value')
#     plt.title(f"Attitude Quaternion: {lab[i]}")
#     plt.xlabel('Time (s)')
#     plt.ylabel(f"{lab[i]}")
#     plt.legend()
# plt.show()

# phi = np.zeros(len(time_arr),)
# th = np.zeros_like(phi)
# psi = np.zeros_like(phi)
# for i in range(len(qu)):
#     phi[i], th[i], psi[i] = quat2euler(qu[:,[i]])

# # q norm
# plt.figure()
# plt.plot(time_arr, qnorm)
# plt.title("Attitude Quaternion Norm")
# plt.xlabel("Time (s)")
# plt.ylabel("norm(q)")
# ax = plt.gca()
# ax.get_yaxis().get_major_formatter().set_useOffset(False)
# plt.show()

# # angular velocity
# for i in range(3):
#     sig = np.sqrt(Pu[3+i, 3+i, :])
#     plt.figure()
#     plt.plot(time_arr, Z_rg[i,:], label='measurement')
#     plt.plot(time_arr, wu[i,:], label='estimate')
#     plt.plot(time_arr, w_sim[i,:], label='true value')
#     upper_bound = wu[i,:] + 2*sig
#     lower_bound = wu[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Angular Velocity: $\omega_{i+1}$")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"$\omega_{i+1}$ (rad/s)")
#     plt.legend()
# plt.show()

# # angular veloctiy error
# for i in range(3):
#     sig = np.sqrt(Pu[3+i, 3+i, :])
#     err = wu[i,:] - w_sim[i,:]
#     zerr = Z_rg[i,:] - w_sim[i,:]
#     plt.figure()
#     plt.plot(time_arr, zerr, label='measurement error')
#     plt.plot(time_arr, err, linestyle='--', label='estimate error')
#     upper_bound = 2*sig
#     lower_bound = -2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Angular velocity Error: $\omega_{i+1}$")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"$\omega_{i+1}$ Error (rad/s)")
#     plt.legend()
# plt.show()

# # position
# for i in range(3):
#     lab = ['X', 'Y', 'Z']
#     sig = np.sqrt(Pu[6+i,6+i,:])
#     plt.figure()
#     plt.plot(time_arr, Z_gps[i,:], label='measurement')
#     plt.plot(time_arr, ru[i,:], label='estimate')
#     plt.plot(time_arr, r_sim[i,:], label='true value')
#     upper_bound = ru[i,:] + 2*sig
#     lower_bound = ru[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Position: {lab[i]}")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"{lab[i]} Coordinate (m)")
#     plt.legend()
# plt.show()

# # position error
# for i in range(3):
#     lab = ['X', 'Y', 'Z']
#     sig = np.sqrt(Pu[6+i,6+i,:])
#     err = ru[i,:] - r_sim[i,:]
#     zerr = Z_gps[i,:] - r_sim[i,:]
#     plt.figure()
#     plt.plot(time_arr, zerr, label='measurement error')
#     plt.plot(time_arr, err, label='estimate error')
#     upper_bound = 2*sig
#     lower_bound = -2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Position Error: {lab[i]}")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"{lab[i]} Coordinate Error (m)")
#     plt.legend()
# plt.show()

# # velocity
# for i in range(3):
#     lab = ['X', 'Y', 'Z']
#     sig = np.sqrt(Pu[9+i,9+i,:])
#     plt.figure()
#     plt.plot(time_arr, vu[i,:], label='estimate')
#     plt.plot(time_arr, v_sim[i,:], label='true value')
#     upper_bound = vu[i,:] + 2*sig
#     lower_bound = vu[i,:] - 2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound, 
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Velocity: {lab[i]}")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"{lab[i]} Coordinate (m/s)")
#     plt.legend()
# plt.show()

# # velocity error
# for i in range(3):
#     lab = ['X', 'Y', 'Z']
#     sig = np.sqrt(Pu[9+i,9+i,:])
#     err = vu[i,:] - v_sim[i,:]
#     plt.figure()
#     plt.plot(time_arr, err, label='estimate error')
#     upper_bound = 2*sig
#     lower_bound = -2*sig
#     plt.fill_between(time_arr, lower_bound, upper_bound,
#                      color='gray', alpha=0.2, label='$\pm 2\sigma$')
#     plt.title(f"Velocity Error: {lab[i]}")
#     plt.xlabel("Time (s)")
#     plt.ylabel(f"{lab[i]} Coordinate Error (m/s)")
#     plt.legend()
# plt.show()
