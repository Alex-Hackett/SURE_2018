# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 13:20:30 2018

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
e_mass = 9.10938356e-31
m_ph = 3 * 10**(-5) * e_mass
m_x = 0.08*e_mass
hbar = 6.582119514e-16

rabi = 26e-3 
omega_0 = 1.7
epsilon_0 = 0
c = 299792458
delta = 5.4e-3


def polariton_dispersion(k_mod):
    sum_term = (delta + ((hbar**2 * k_mod**2)/(2*m_x)) + ((hbar**2 * k_mod**2)/(2*m_ph)))
    diff_term = (delta + ((hbar**2 * k_mod**2)/(2*m_x)) - ((hbar**2 * k_mod**2)/(2*m_ph)))
    e_cav = ((hbar**2 * k_mod**2)/(2*m_ph)) + delta
    E_l = 0.5 * ((sum_term) - (np.sqrt(diff_term**2 + rabi**2)))
    e_ex = ((hbar**2 * k_mod**2)/(2*m_x))
    E_p = 0.5 * ((sum_term) + (np.sqrt(diff_term**2 + rabi**2)))
    return E_l, E_p, e_cav, e_ex



ks = np.linspace(10**(-3), 10**(-1), 10000)
ks = ks * (1e9)
energysl, energysp, cav, ex = polariton_dispersion(ks)
ks = ks * (1e-9)
fig1 = plt.figure()
plt.semilogx(ks, energysl*1e3, label = 'Lower Polariton Branch')
plt.semilogx(ks, energysp*1e3, label = 'Upper Polariton Branch')
plt.semilogx(ks, cav*1e3, label = 'Photon')
plt.semilogx(ks, ex*1e3, label = 'Exciton')
plt.title('Polariton Dispersion Curves')
plt.xlabel(r'Wavevector [$10^7 cm^{-1}$]')
plt.ylabel(r'Energy [meV]')
plt.legend()
plt.ylim(-15,40)
plt.xlim(10**(-3))
plt.show()

    
    
    
    