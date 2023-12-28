import user_parameters.user_inputs  as user

import cmath
import math
import numpy as np
import pandas as pd

def import_data(csv_file_name):
    """
    Import data from a CSV file and extract relevant parameters.
    
    Parameters:
    - csv_file_name (str): The name of the CSV file.
    
    Returns:
    - alpha_pn_dB (numpy.ndarray): Array of PN-junction loss values in dB/cm.
    - neff_pn (numpy.ndarray): Array of PN-junction effective index values.
    - V (numpy.ndarray): Array of V_anode values, assuming Vcathode = 0.
    """
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file_name)

    # Extract relevant parameters
    alpha_pn_dB = np.squeeze(df['loss_TE']) * 1e-2
    neff_pn = np.squeeze(df['neff_TE'])
    V = np.squeeze(df['V_anode'])

    return alpha_pn_dB, neff_pn, V

def calculate_round_trip_loss_phase(alpha_i_dB, alpha_pn_dB, neff_i, neff_pn, wavelength):
    """
    Calculate total phase accumulation, total optical attenuation, and total length.
    
    Parameters:
    - alpha_i_dB (float): Intrinsic loss in dB/cm.
    - alpha_pn_dB (float): PN-junction loss in dB/cm.
    - neff_i (float): Intrinsic effective index.
    - neff_pn (float): PN-junction effective index.
    - wavelength (float): Wavelength in meters.
    - Lc (float): Length of undoped section in meters.
    - R (float): Radius of the ring in meters.
    - L_pn (float): Length of the PN-junction in meters.
    
    Returns:
    - phi_rt (float): Total phase accumulation.
    - A (float): Total optical attenuation.
    - L_rt (float): Total length.
    """

    # Converts the intrinsic and pn-junction loss from dB/cm to 1/cm
    alpha_i = -math.log(10**(-alpha_i_dB / 10))
    alpha_pn = -math.log(10**(-alpha_pn_dB / 10))

    # Separates the Intrinsic and PN-junction loss for any undoped sections
    alpha_pn += alpha_i
    alpha_exc = alpha_i 

    # Calculates the doped RT and un-doped lengths
    L_rt = 2 * user.Lc + 2 * user.PI * user.R
    L_exc = L_rt - user.L_pn

    # Calculates the phases due to doped RT and un-doped lengths    
    phi_pn = (2 * user.PI / wavelength) * neff_pn * user.L_pn
    phi_exc = (2 * user.PI / wavelength) * user.neff_i * L_exc
    
    # Total Phase Accumulation
    phi_rt = phi_pn + phi_exc

    # Calculates the loss due to doped RT and un-doped lengths    
    A_pn = math.exp(-alpha_pn * 100 * user.L_pn)  
    A_exc = math.exp(-alpha_exc * 100 * L_exc)                     

    # Total Optical Attenuation and RT total loss
    A = A_pn * A_exc 
    A_rt = np.exp(-A*L_rt)

    return phi_rt, A, L_rt, A_rt


def calculate_electric_fields(k1, k2, A, L_rt, phi_rt, wavelength, filter_type, Ng):
    """
    Calculate optical fields (Ethru and Edrop) and coupled Q-factor.
    
    Parameters:
    - k1 (float): Coupling Coefficient for Ethru.
    - k2 (float): Coupling Coefficient for Edrop.
    - A (float): Total optical attenuation.
    - L_rt (float): Total length in meters.
    - phi_rt (float): Total phase accumulation in rad.
    - wavelength (float): Wavelength in meters.
    - filter_type (str): Type of optical filter ('all-pass' or 'add-drop').
    - Ng (float): Group index.
    
    Returns:
    - Ethru (complex): Optical field for the through port.
    - Edrop (complex): Optical field for the drop port.
    - Qc (float): Coupled Q-factor.
    """
    
    if filter_type == 'all-pass':
        t = cmath.sqrt(1 - k1**2)
        Ethru = (-cmath.sqrt(A) + t * cmath.exp(-1j * phi_rt)) / (-cmath.sqrt(A) * np.conj(t) + cmath.exp(-1j * phi_rt))
        Edrop = 0
        Qc = -(math.pi * L_rt * Ng) / (wavelength * math.log(abs(t)))
    elif filter_type == 'add-drop':
        t1 = cmath.sqrt(1 - k1**2)
        t2 = cmath.sqrt(1 - k2**2)
        Ethru = (t1 - np.conj(t2) * cmath.sqrt(A) * cmath.exp(1j * phi_rt)) / (
                    1 - cmath.sqrt(A) * np.conj(t1) * np.conj(t2) * cmath.exp(1j * phi_rt))
        Edrop = -np.conj(k1) * k2 * cmath.sqrt(cmath.sqrt(A)) * cmath.exp(1j * phi_rt / 2) / (
                    1 - cmath.sqrt(A) * np.conj(t1) * np.conj(t2) * cmath.exp(1j * phi_rt))
        Qc1 = -(math.pi * L_rt * Ng) / (wavelength * math.log(abs(t1)))
        Qc2 = -(math.pi * L_rt * Ng) / (wavelength * math.log(abs(t2)))
        Qc = 1 / (1 / Qc1 + 1 / Qc2)
    else:
        raise ValueError("The 'filter_type' has to be 'all-pass' or 'add-drop'.\n")
        
    return Ethru, Edrop, Qc

def calculate_radius_at_resonance(wavelength, neff_i, L_rt):
    """
    Calculate ring radius for resonance given wavelength, effective index, and total length.

    Parameters:
    - wavelength (float): Resonance wavelength in meters.
    - neff_i (float): Intrinsic effective index.
    - L_rt (float): Total length in meters.

    Returns:
    - R (float): Ring radius in meters.
    """    
    m = math.floor(neff_i * L_rt / wavelength)
    L_rt_new = wavelength * m / neff_i
    R = L_rt_new / (2 * user.PI)
    wavelength_new = user.neff_i * L_rt / m

    return R, wavelength_new

def convert_FSR_to_radius(FSR, Ng):
    """
    Convert Free Spectral Range (FSR) to ring radius.
    
    Parameters:
    - FSR (float): Free Spectral Range in Hertz.
    - Ng (float): Group index.
    - c (float): Speed of light in meters per second. Default is 3e8 m/s.
    
    Returns:
    - R (float): Ring radius in meters.
    """
    L_rt = user.c / (Ng*FSR)
    R = L_rt / (2*user.PI)

    return R


def critical_coupling(k2, A_rt, filter_type):
    """
    Calculate critical coupling for different filter types.

    Parameters:
    - k2 (float): Coupling coefficient for the filter.
    - A_rt (float): Total Round Trip loss .
    - filter_type (str): Type of filter, either 'all-pass' or 'add-drop'.

    Returns:
    - k1 (float): Calculated critical coupling coefficient.
    """
    # Default value for k1
    k1 = None
    
    if filter_type == 'all-pass':
        # No action needed for 'all-pass' filter type
        print("No action needed for 'all-pass' filter type.")
    elif filter_type == 'add-drop':
        # Calculate k1 for 'add-drop' filter type
        k1 = A_rt * k2
    else:
        # Handle unexpected filter types gracefully
        print(f"Warning: Unexpected filter_type '{filter_type}'. Defaulting to None.")
    
    return k1

if __name__ == "__main__":
    # Calculate Round Trip Loss Phase
    A, phi_rt, L_rt, A_rt = calculate_round_trip_loss_phase(
        user.alpha_i_dB, user.alpha_pn_dB, user.neff_i, user.neff_pn, user.wavelength_0)
    print(f'Loss [1/cm]: {round(A, 3)}\n'
          f'Total Round Trip Phase [rad]: {round(phi_rt, 3)}\n'
          f'Round Trip Length [um]: {round(L_rt * 1e6, 3)}\n'
          f'Total Round Trip Loss: {round(A_rt,5)}\n')

    # Calculate Electric Fields
    Ethru, Edrop, Qc = calculate_electric_fields(
        user.k1, user.k2, A, L_rt, phi_rt, user.wavelength_0, user.filter_type, user.Ng)
    print(f'Coupled Q-factor: {round(Qc)}\n')

    # Calculate Radius at Resonance
    R, wavelength_new = calculate_radius_at_resonance(
        user.wavelength_0, user.neff_i, L_rt)
    print(f'Calculate radius at resonance wavelength {user.wavelength_0 * 1e6} um\n'
          f'Radius [um]: {round(R * 1e6, 3)}\n')

    # Convert FSR to Radius
    R = convert_FSR_to_radius(FSR=user.FSR, Ng=user.Ng)
    print(f'Calculate radius with FSR of {round(user.FSR * 1e-12, 3)} THz\n'
          f'Radius [um]: {round(R * 1e6, 3)}\n')

    