# -*- coding: utf-8 -*-
"""
Created on Fri May 10 13:50:38 2019

@author: Harald
"""
import numpy as np

# Mixing func:
def NB_Mix(HPAA, HPAa, HPtot):
    
    if HPtot == 0:
        HPtot = 1e-14
    
    mix = (HPAA + HPAa/2)/HPtot
    return mix

# Host func:
def NB_H(Ho, Hotot, Po, lam, K, a):
    
    if Ho < 1e-14:
        H = 0
    else:
        H = Ho*lam*np.exp(-Hotot*np.log(lam)/K)*np.exp(-a*Po)
    
    return H

# Parasite func:
def NB_P(Ho, Po, a, c):
    
    if Po < 1e-14:
        P = 0
    else:
        P = c*Ho*(1 - np.exp(-a*Po))
    
    return P