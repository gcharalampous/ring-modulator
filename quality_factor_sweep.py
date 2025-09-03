from ring_utils import critical_coupling, calculate_round_trip_loss_phase, calculate_electric_fields, load_user_config, load_user_config
import numpy as np
import matplotlib.pyplot as plt



if __name__ == "__main__":
   # Calculate Round Trip Loss Phase
    user = load_user_config()
    phi_rt, A, L_rt, A_rt = calculate_round_trip_loss_phase(
                alpha_i_dB_per_cm=user.alpha_i_dB_per_cm, 
                alpha_pn_dB_per_cm=user.alpha_pn_dB_per_cm,
                neff_i=user.neff_i, 
                neff_pn=user.neff_pn, 
                wavelength=user.wavelength_0,
                L_pn=user.L_pn, Lc=user.Lc, radius=user.R)
    
    # Vectorize coupling coefficients
    Qc = np.zeros(user.pts)
    Ql = np.zeros(user.pts)
    k2 = np.linspace(0.01,0.99,user.pts)
    k1 = np.zeros(user.pts)
    px = 1 / plt.rcParams['figure.dpi']
    
    for i in range(0, user.pts):
        k1[i] = critical_coupling(k2 = k2[i], A_rt=A_rt, filter_type = user.filter_type)
        Ethru, Edrop, Qc[i], Ql[i] = calculate_electric_fields(
        k1[i], k2[i], A, L_rt, phi_rt, user.wavelength_0, user.filter_type, user.Ng)
    
    # Plotting |k1|^2 vs |k2|^2
    plt.figure(1, figsize=(512 * px, 256 * px))
    plt.loglog(k1**2,k2**2)
    plt.xlabel('$|k_1|^2$')
    plt.ylabel("$|k_2|^2$")
    plt.grid(which='both')
    plt.tight_layout()

    # Plotting |k1|^2 vs Qc
    plt.figure(2, figsize=(512 * px, 256 * px))
    plt.loglog(k1**2,Qc)
    plt.ylabel('Coupled Q-factor')
    plt.xlabel("$|k_1|^2$")
    plt.grid(which='both')
    plt.tight_layout()

    # Show the plots
    plt.show()