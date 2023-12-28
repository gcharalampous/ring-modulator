import user_parameters.user_inputs  as user

import matplotlib.pyplot as plt
import numpy as np
from functions import import_data, calculate_round_trip_loss_phase, calculate_electric_fields

# Input data
alpha_pn_dB_, neff_pn_, V = import_data(user.csv_voltage_neff_sweep)
N = len(alpha_pn_dB_)
wavelength_array = np.linspace(user.wavelength_0 - user.wavelength_bw/2, 
                               user.wavelength_0 + user.wavelength_bw/2, user.pts)

# Arrays for storing results
phi_rt = [0] * len(wavelength_array)
A = [0] * len(wavelength_array)
Qc = [0] * len(wavelength_array)
Edrop_max = [0] * N
Ethru_min = [0] * N
wave_max = [0] * N

# Initialize 2D arrays (wavelength, Voltage)
rows, cols = len(wavelength_array), N
Edrop = np.zeros((rows, cols), dtype=complex)
Ethru = np.zeros((rows, cols), dtype=complex)
Qc = np.zeros((rows, cols), dtype=float)

if __name__ == "__main__":
    # Get the pixel size for figure dimensions
    px = 1 / plt.rcParams['figure.dpi']

    # Loop through each voltage value
    for j, voltage in enumerate(V):
        # Loop through each wavelength value
        for i, wavelength in enumerate(wavelength_array):
            # Calculate loss and phase
            phi_rt, A, L_rt, A_rt = calculate_round_trip_loss_phase(
                alpha_i_dB=alpha_pn_dB_[0], alpha_pn_dB=alpha_pn_dB_[j],
                neff_i=np.real(complex(neff_pn_[0])), neff_pn=np.real(complex(neff_pn_[j])), wavelength = wavelength)

            # Calculate optical fields
            Ethru[i, j], Edrop[i, j], Qc[i, j] = calculate_electric_fields(user.k1, user.k2, A, L_rt, phi_rt, wavelength, filter_type=user.filter_type, Ng=user.Ng)

        # Record max drop and min through values, map to wavelength (use the data to calculate the tunability rate)
        Edrop_max[j] = np.max(20 * np.log10(abs(Edrop[:, j])))
        Ethru_min[j] = np.min(20 * np.log10(abs(Ethru[:, j])))
        wave_max[j] = wavelength_array[np.argmax(20 * np.log10(abs(Edrop[:, j])))]

        # Plot Through-port
        plt.figure(1, figsize=(780 * px, 256 * px))
        plt.plot(wavelength_array * 1e9, 20 * np.log10(abs(Ethru[:, j])), label=f"{voltage} V")
        plt.figure(3, figsize=(780 * px, 256 * px))
        plt.plot(wavelength_array * 1e9, np.angle(Ethru[:, j]), label=f"{voltage} V")

        # Plot Add-port if filter_type is 'add-drop'
        if user.filter_type == 'add-drop':
            plt.figure(2, figsize=(780 * px, 256 * px))
            plt.plot(wavelength_array * 1e9, 20 * np.log10(abs(Edrop[:, j])), label=f"{voltage} V")
            plt.figure(4, figsize=(780 * px, 256 * px))
            plt.plot(wavelength_array * 1e9, np.angle(Edrop[:, j]), label=f"{voltage} V")
    # Common settings for both plots
    for i in range(1, 3):
        plt.figure(i)
        plt.title('Through-port' if i == 1 else 'Add-port')
        plt.grid()
        plt.legend()
        plt.ylabel("T (dB)")
        plt.xlabel("Wavelength (nm)")
    for i in range(3, 5):
        plt.figure(i)
        plt.title('Through-port' if i == 3 else 'Add-port')
        plt.grid()
        plt.legend()
        plt.ylabel("Angle (rad)")
        plt.xlabel("Wavelength (nm)")

    # Show the plots
    plt.show()
