# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:29:04 2018

@author: Alex
"""

import os
import astropy as ap
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import math
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import sys
import time
from random import randrange, random, choice
import os.path



def LP_dis(mod_k):
    hbar = 6.582119514e-16 #in eV.s
    lz = 300e-6 #Cavity Lenght in Meters
    M = 1 #Mode of coupled Light
    c = 3e8
    n0 = 3.9476 #Refective index of CAVITY
    qz = (2* M * np.pi) / (lz)
    E_ex_0 = 1.557
    rabi = 30 * 1e-3#Rabi splitting in eV
    
    
    E_cav = hbar * c/n0 * np.sqrt(qz**2 + mod_k**2) #Cavity Mode dispersion
    E_ex = E_ex_0 #Exciton Dispersion
    d_k = E_cav - E_ex
    #Lower Branch Polariton Dispersion
    E_l = (E_cav + E_ex)/2 + (np.sqrt(d_k**2 + 4 * rabi**2))/2
    return E_l

k_array = np.arange(11) * 1e6
omega_array = np.zeros(11)
j = 0
for i in k_array:
    omega_array[j] = LP_dis(i)
    j += 1

plt.plot(k_array*1e-6,omega_array * 1e3)
    



