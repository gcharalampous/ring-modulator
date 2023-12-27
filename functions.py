#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# @author: Georgios Gcharalampous (gcharalampous)
# version ='1.0'
# ---------------------------------------------------------------------------


from user_parameters.user_inputs import *


import cmath
import math
import numpy as np
import pandas as pd

def import_data(csv_file_name):
    df = pd.read_csv(csv_file_name)

    alpha_pn_dB = np.squeeze(df['loss_TE'])*1e-2
    neff_pn = np.squeeze(df['neff_TE'])
    V = np.squeeze(df['V_anode'])

    return alpha_pn_dB, neff_pn, V


def calculate_loss_RTphase(alpha_i_dB, alpha_pn_dB, neff_i, neff_pn, wavelength):
    
    # Converts the intrinsic and pn-junction loss from dB/cm to 1/cm
    alpha_i = -math.log(10**(-alpha_i_dB/10))
    alpha_pn = -math.log(10**(-alpha_pn_dB/10))


    # Separates the Intrinsic and PN-junction loss for any undoped sections
    alpha_pn = alpha_i + alpha_pn
    alpha_exc = alpha_i 


    # Calculates the doped RT and un-doped lengths
    L_rt = 2*Lc + 2*PI*R
    L_exc = L_rt - L_pn

    # Calculates the phases due to doped RT and un-doped lengths    
    phi_pn = (2*PI/wavelength)*neff_pn*L_pn
    phi_exc = (2*PI/wavelength)*neff_i*L_exc
    
    # Total Phase Accumulation
    phi_rt = phi_pn + phi_exc

    # Calculates the loss due to doped RT and un-doped lengths    
    A_pn = math.exp(-alpha_pn*100*L_pn);        
    A_exc = math.exp(-alpha_exc*100*L_exc);                     

    # Total Optical Attenuation
    A = A_pn * A_exc; 

    return phi_rt, A, L_rt




def ring_optical_fields(k1,k2,A, L_rt, phi_rt, wavelength):
    if filter_type == 'all-pass':
        t = cmath.sqrt(1 - k1**2)
        Ethru = (-cmath.sqrt(A) + t * cmath.exp(-1j * phi_rt)) / (-cmath.sqrt(A) * np.conj(t) + cmath.exp(-1j * phi_rt))
        Edrop = 0
        Qc = -(math.pi * L_rt * Ng) / (wavelength * math.log(abs(t)))

        # print(filter_type)
        # print("Coupled Q-factor: " + str(Qc))


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

        # print(filter_type)
        # print("Coupled Q-factor: " + str(Qc))

    else:
        raise ValueError("The 'filter_type' has to be 'all-pass' or 'add-drop'.\n")
    
    return Ethru, Edrop, Qc



def resonance_wavelength_to_Radii(wavelength, neff_i, L_rt):
    m = neff_i*L_rt/wavelength
    m = math.floor(m)
    L_rt_new = wavelength*m/neff_i
    R = L_rt_new/(2*PI)

    wavelength_new = neff_i*L_rt/m
    print('R='+str(R) + 'at Resonance '+str(wavelength_new*1e6)+'um' )
    return R

def FSR_to_Radii(FSR, Ng):
    L_rt = c/FSR
    R = L_rt/(2*PI)
    print('R='+str(R) + 'at '+str(FSR*1e-12)+'THz' )
    return R

if(__name__ == "__main__"):

    A, phi_rt, L_rt = calculate_loss_RTphase(alpha_i_dB, alpha_pn_dB, neff_i, neff_pn, wavelength)
    
    Ethru, Edrop, Qc = ring_optical_fields(k1,k2,A, L_rt, phi_rt, wavelength)

    R = resonance_wavelength_to_Radii(wavelength, neff_i, L_rt)