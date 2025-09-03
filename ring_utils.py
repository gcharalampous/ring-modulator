import yaml
import pandas as pd
from pathlib import Path
from types import SimpleNamespace
import math
import cmath
import numpy as np
import scipy.constants as const

def load_user_config():
    """
    Load user parameters from YAML and return as SimpleNamespace.
    """


    BASE_DIR = Path(__file__).resolve().parent
    CONFIG_PATH = BASE_DIR / 'user_parameters' / 'user_inputs.yaml'
    with open(CONFIG_PATH, 'r') as f:
        cfg = yaml.safe_load(f)

    # Resolve CSV path relative to the package root (BASE_DIR)
    csv_path = cfg.get('csv_voltage_neff_sweep')
    if csv_path is not None:
        csv_p = (BASE_DIR / csv_path).resolve()
        cfg['csv_voltage_neff_sweep'] = str(csv_p)

    # Evaluate simple derived expressions if present (e.g., L_pn: "2*PI*R")
    try:
        if isinstance(cfg.get('L_pn'), str):
            locals_map = {k: v for k, v in cfg.items() if isinstance(v, (int, float))}
            if 'PI' in cfg:
                locals_map['PI'] = cfg['PI']
            if 'R' in cfg:
                locals_map['R'] = cfg['R']
            cfg['L_pn'] = float(eval(cfg['L_pn'], {}, locals_map))
    except Exception:
        pass

    return SimpleNamespace(**cfg)

# Usage: user = load_user_config()

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
    df = pd.read_csv(csv_file_name, comment="#")

    # Extract relevant parameters
    alpha_pn_dB_per_cm = np.squeeze(df['loss_TE']) * 1e-2
    neff_pn = np.squeeze(df['neff_TE'])
    V = np.squeeze(df['V_anode'])

    return alpha_pn_dB_per_cm, neff_pn, V

def calculate_round_trip_loss_phase(alpha_i_dB_per_cm, alpha_pn_dB_per_cm, L_pn, Lc, radius, neff_i, neff_pn, wavelength):
    """
    Calculate total phase accumulation, total optical attenuation, and total length.

    Parameters:
    - alpha_i_dB_per_cm (float): Intrinsic loss in dB/cm.
    - alpha_pn_dB_per_cm (float): PN-junction loss in dB/cm.
    - L_pn (float): Length of the PN-junction in meters.
    - Lc (float): Length of undoped section in meters.
    - radius (float): Radius of the ring in meters.
    - neff_i (float): Intrinsic effective index.
    - neff_pn (float): PN-junction effective index.
    - wavelength (float): Wavelength in meters.

    Returns:
    - phi_rt (float): Total phase accumulation.
    - A (float): Total optical attenuation.
    - L_rt (float): Total round-trip length.
    - A_rt (float): Total round-trip loss.
    """

    # Converts the intrinsic and pn-junction loss from dB/cm to 1/cm
    alpha_i_per_cm = -math.log(10**(-alpha_i_dB_per_cm / 10))
    alpha_pn_per_cm = -math.log(10**(-alpha_pn_dB_per_cm / 10))

    # Separates the Intrinsic and PN-junction loss for any undoped sections
    alpha_pn_per_cm += alpha_i_per_cm
    alpha_exc = alpha_i_per_cm 

    # Calculates the doped RT and un-doped lengths
    L_rt = 2 * Lc + 2 * math.pi * radius
    L_exc = L_rt - L_pn

    # Calculates the phases due to doped RT and un-doped lengths    
    phi_pn = (2 * math.pi / wavelength) * neff_pn * L_pn
    phi_exc = (2 * math.pi / wavelength) * neff_i * L_exc
    
    # Total Phase Accumulation
    phi_rt = phi_pn + phi_exc

    # Calculates the loss due to doped RT and un-doped lengths    
    A_pn = math.exp(-alpha_pn_per_cm * 100 * L_pn)  
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
    - Qloaded (float): Loaded Q-factor.
    """
    # Find notch
    if filter_type == 'all-pass':
        t = cmath.sqrt(1 - k1**2)
        Ethru = (-cmath.sqrt(A) + t * cmath.exp(-1j * phi_rt)) / (-cmath.sqrt(A) * np.conj(t) + cmath.exp(-1j * phi_rt))
        Edrop = 0
        Qc = -(math.pi * L_rt * Ng) / (wavelength * math.log(abs(t)))
        Qloaded = (np.pi * Ng * L_rt * np.sqrt(t * A)) / (wavelength * (1 - t * A))
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
        Qloaded = (np.pi * Ng * L_rt * np.sqrt(t1 * t2 * A)) / (wavelength * (1 - t1 * t2 * A))
    else:
        raise ValueError("The 'filter_type' has to be 'all-pass' or 'add-drop'.\n")
       
    return Ethru, Edrop, Qc, Qloaded

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
    m_raw = neff_i * L_rt / wavelength
    m = max(1, round(m_raw))  # Ensure m is at least 1
    L_rt_new = wavelength * m / neff_i
    R = L_rt_new / (2 * math.pi)
    wavelength_new = neff_i * L_rt / m

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
    L_rt = const.c / (Ng*FSR)
    R = L_rt / (2*const.pi)

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
        k1 = math.sqrt(1-((A_rt**2) * (1-k2**2)))
    else:
        # Handle unexpected filter types gracefully
        print(f"Warning: Unexpected filter_type '{filter_type}'. Defaulting to None.")
    
    return k1

def extinction_ratio_vs_wavelength(E_num, E_den, wavelength, mode="power", floor=1e-18):
    """
    Compute extinction ratio ER(λ) between two complex fields across wavelength.

    Parameters
    ----------
    E_num : array-like of complex
        Numerator field vs wavelength (e.g., E_through or E_drop).
    E_den : array-like of complex
        Denominator field vs wavelength (e.g., reference/baseline).
    wavelength : array-like of float
        Wavelength samples (same length as E_num/E_den).
    mode : {"power", "field"}
        "power": ER_dB = 10*log10(|E_num|^2 / |E_den|^2)  [recommended]
        "field": ER_dB = 20*log10(|E_num| / |E_den|)
    floor : float
        Magnitude floor to avoid divide-by-zero / log of zero.

    Returns
    -------
    er_db : np.ndarray
        Extinction ratio in dB at each wavelength.
    er_lin : np.ndarray
        Extinction ratio in linear scale at each wavelength.
    H : np.ndarray of complex
        Complex transfer ratio H = E_num / E_den (useful for phase, etc.).
    """
    E_num = np.asarray(E_num)
    E_den = np.asarray(E_den)
    wl = np.asarray(wavelength, float)

    if E_num.shape != E_den.shape or E_num.shape != wl.shape:
        raise ValueError("E_num, E_den, and wavelength must have the same 1D shape.")

    # Complex ratio (keeps amplitude & phase info)
    # Add tiny floor to denominator magnitude to avoid numerical blow-ups
    denom = E_den.copy()
    small = (np.abs(denom) < floor)
    denom[small] = denom[small] + floor  # nudge tiny/zero values

    H = E_num / denom
    mag_num = np.maximum(np.abs(E_num), floor)
    mag_den = np.maximum(np.abs(E_den), floor)

    if mode == "power":
        er_lin = (mag_num**2) / (mag_den**2)
        er_db = 10.0 * np.log10(er_lin)
    elif mode == "field":
        er_lin = mag_num / mag_den
        er_db = 20.0 * np.log10(er_lin)
    else:
        raise ValueError("mode must be 'power' or 'field'.")

    return er_db.astype(float), er_lin.astype(float), H



if __name__ == "__main__":
    # Load user configuration
    user = load_user_config()

    # Calculate Round Trip Loss Phase
    phi_rt, A, L_rt, A_rt = calculate_round_trip_loss_phase(
        user.alpha_i_dB_per_cm, user.alpha_pn_dB_per_cm, user.L_pn, user.Lc, user.R, user.neff_i, user.neff_pn, user.wavelength_0)
    print(f'Loss [1/cm]: {round(A, 3)}\n'
          f'Total Round Trip Phase [rad]: {round(phi_rt, 3)}\n'
          f'Round Trip Length [um]: {round(L_rt * 1e6, 3)}\n'
          f'Total Round Trip Loss: {round(A_rt,5)}\n')

    # Calculate Critical Coupling
    k1 = critical_coupling(user.k2, A_rt, user.filter_type)
    if k1 is not None:
        print(f'Critical Coupling [1]: {round(k1, 3)}\n')
    else:
        print('Critical Coupling [1]: None\n')

    # Calculate Electric Fields
    Ethru, Edrop, Qc, Qloaded = calculate_electric_fields(
        user.k1, user.k2, A, L_rt, phi_rt, user.wavelength_0, user.filter_type, user.Ng)
    print(f'Coupled Q-factor: {round(abs(Qc), 3)}\n'
          f'Loaded Q-factor: {round(abs(Qloaded), 3)}\n')

    # Calculate Radius at Resonance
    R, wavelength_new = calculate_radius_at_resonance(
        user.wavelength_0, user.neff_i, L_rt)
    print(f'Calculate radius at resonance wavelength {user.wavelength_0 * 1e6} um\n'
          f'Radius [um]: {round(R * 1e6, 3)}\n')

    # Convert FSR to Radius
    R = convert_FSR_to_radius(FSR=user.FSR, Ng=user.Ng)
    print(f'Calculate radius with FSR of {round(user.FSR * 1e-12, 3)} THz\n'
          f'Radius [um]: {round(R * 1e6, 3)}\n')

    
