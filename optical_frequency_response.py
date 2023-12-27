
from user_parameters.user_inputs import *
from functions import *
import matplotlib.pyplot as plt
import numpy as np

csv_file_name = 'input_data/neff_Voltage_sweep_TE.csv'


alpha_pn_dB_, neff_pn_, V = import_data(csv_file_name)
N = len(alpha_pn_dB_)
wavelength_array = np.linspace(1.54e-6,1.56e-6,500)
phi_rt=[0]*len(wavelength_array)
A=[0]*len(wavelength_array)
Qc = [0]*len(wavelength_array)
Edrop_max = [0]*N
Ethru_min = [0]*N
wave_max = [0]*N


# Initialize 2D arrays (wavelength,Voltage)
rows = len(wavelength_array)
cols = N
Edrop = np.zeros((rows, cols), dtype=complex)
Ethru = np.zeros((rows, cols), dtype=complex)
Qc = np.zeros((rows, cols), dtype=float)


if(__name__ == "__main__"):



    px = 1/plt.rcParams['figure.dpi']  # pixel in inches

    for j in range(N):
        i=0
        for w_element in wavelength_array:
            

            phi_rt, A, L_rt = calculate_loss_RTphase(
                alpha_i_dB = alpha_pn_dB_[0], alpha_pn_dB = alpha_pn_dB_[j],
                neff_i = np.real(complex(neff_pn_[0])), neff_pn = np.real(complex(neff_pn_[j])), wavelength=w_element)

            Ethru[i,j], Edrop[i,j], Qc[i,j] = ring_optical_fields(k1,k2,A, L_rt, phi_rt, wavelength = w_element)
            i = i + 1

        Edrop_max[j] = np.max(20*np.log10(abs(Edrop[:,j])))
        Ethru_min[j] = np.min(20*np.log10(abs(Ethru[:,j])))
        wave_max[j] = wavelength_array[np.argmax(20*np.log10(abs(Edrop[:,j])))]
        plt.figure(1,figsize=(780*px, 256*px))
        plt.plot(wavelength_array*1e9,20*np.log10(abs(Ethru[:,j])), label = str(V[j]) + 'V')
        plt.title('Through-port')
        plt.grid()
        plt.legend()
        plt.ylabel("T (dB)")
        plt.xlabel("wavelength (nm)")
        
        if(filter_type == 'add-drop'):
    
            plt.figure(2,figsize=(780*px, 256*px))
            plt.plot(wavelength_array*1e9,20*np.log10(abs(Edrop[:,j])), label = str(V[j])+' V')
            plt.title('Add-port')
            plt.grid()
            plt.legend()
            plt.ylabel("T (dB)")
            plt.xlabel("wavelength (nm)")