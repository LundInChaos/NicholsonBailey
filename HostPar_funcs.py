# -*- coding: utf-8 -*-
"""
Created on Fri May 10 13:50:38 2019

@author: Harald
"""
import numpy as np

# Mixing func:
# A function that is used to calculate how the genotypes of host/parasite
# gets distributed after mixing. When calculating for hosts only host
# populations are used, for parasites only parasite populations.
def NB_Mix(HPAA, HPAa, HPtot):
    
    # Avoid division by zero later:
    if HPtot == 0:
        HPtot = 1e-14
    
    mix = (HPAA + HPAa/2)/HPtot
    return mix

# Host func:
# A function that calculates the host population (before mixing)
# version considering genotypes of hosts
def NB_H(Ho, Hotot, Po, lam, K, a):
    
    # Pop goes extinct if density is too low:
    if Ho < 1e-14:
        H = 0
    else:
        H = Ho*lam*np.exp(-Hotot*np.log(lam)/K)*np.exp(-a*Po)
    
    return H

# Parasite func:
# A function that calculates the parasite population (before mixing)
def NB_P(Ho, Po, a, c):
    
    # Pop goes extinct if density is too low:
    if Po < 1e-14:
        P = 0
    else:
        P = c*Ho*(1 - np.exp(-a*Po))
    
    return P