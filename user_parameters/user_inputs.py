# CONSTANTS
c = 299792458                           # Speed of light [m/s]
PI = 3.141592653589793

# ---------------------------------------

# External files to input
csv_voltage_neff_sweep = 'input_data/neff_Voltage_sweep_TE.csv'


# SIMULATION PARAMETERS
wavelength_0 = 1.55e-6                  # Single wavelength point
wavelength_bw = 0.02e-6                 # Wavelength Bandwidth
pts = 512                               # Number of points (wavelength)

# ---------------------------------------


# RING MODULATOR PARAMETERS
filter_type = 'add-drop'                # Change this to 'all-pass' for the other case
k1 = 0.18                               # Coupling Coefficient k1
k2 = 0.23                               # Coupling Coefficient k2
R = 3.95e-6                             # Ring Resonator Radius
FSR = 3.2e12                            # Desires FSR in Hz
Lc = 0                                  # Straight Coupling Length
L_pn = 2*PI*R                           # Total Doped Ring circumferance

# ---------------------------------------


# WAVEGUIDE PARAMETERS
Ng = 4.029
neff_i = 2.4998                         # Effective Index of intrinsic region at wavelength_0
neff_pn = 2.4969                        # Effective Index of PN-Junction at wavelength_0
alpha_i_dB = 5                          # Intrinsic waveguide optical loss, in dB/cm
alpha_pn_dB = 70                        # PN-junction waveguide optical loss, in dB/cm   

# ---------------------------------------




















