from functions import critical_coupling, calculate_round_trip_loss_phase, calculate_electric_fields
import user_parameters.user_inputs  as user

import numpy as np
import matplotlib.pyplot as plt



if __name__ == "__main__":
   # Calculate Round Trip Loss Phase
    A, phi_rt, L_rt, A_rt = calculate_round_trip_loss_phase(
        user.alpha_i_dB, user.alpha_pn_dB, user.neff_i, user.neff_pn, user.wavelength_0)
    
    # Vectorize coupling coefficients
    Qc = np.zeros(user.pts)
    k2 = np.linspace(0.01,0.99,user.pts)
    k1 = np.zeros(user.pts)
    px = 1 / plt.rcParams['figure.dpi']
    
    for i in range(0, user.pts):
        k1[i] = critical_coupling(k2 = k2[i], A_rt=A_rt, filter_type = user.filter_type)
        Ethru, Edrop, Qc[i] = calculate_electric_fields(
        k1[i], k2[i], A, L_rt, phi_rt, user.wavelength_0, user.filter_type, user.Ng)
    
    # Plotting |k1|^2 vs |k2|^2
    plt.figure(1, figsize=(512 * px, 256 * px))
    plt.loglog(k1**2,k2**2)
    plt.xlabel('$|k_1|^2$'), plt.ylabel("$|k_2|^2$"), plt.grid(which='both')

    # Plotting |k1|^2 vs Qc
    plt.figure(2, figsize=(512 * px, 256 * px))
    plt.loglog(k1**2,Qc)
    plt.ylabel('Coupled Q-factor'), plt.xlabel("$|k_1|^2$"), plt.grid(which='both')

    # Show the plots
    plt.show()