# -*- coding: utf-8 -*-

import scipy as sp
import os
import astropy as ap
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import math
import numpy as np
import scipy.integrate as integrate
import scipy.misc as spmisc
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import sys
import time
from random import randrange, random, choice
import os.path
import emailer
import polariton_scattering_package as psp
#---------------------------------------------------------
temp = 50
k = 0
k_prime_x = np.linspace(-1e6,1e6,5000)
k_prime_y = np.linspace(-1e6,1e6,5000)
rates = np.zeros((len(k_prime_x), len(k_prime_y)), dtype = float)

for i in range(len(k_prime_x)):
    for j in range(len(k_prime_y)):
        rates[i, j] = psp.W_k_k_prime(k, np.sqrt((k_prime_x[i])**2 + (k_prime_y[j])**2), temp)
    print(i,'/',5000)

rates = np.nan_to_num(rates)

fig1 = plt.figure()
plt.pcolor(k_prime_x*1e-6, k_prime_y*1e-6, rates)
cbar = plt.colorbar()
cbar.set_label(r'Scattering Rate ($s^{-1}$)')
plt.xlabel(r'$k_{x}$ ($\mu m^{-1}$)')
plt.ylabel(r'$k_{y}$ ($\mu m^{-1}$)')
plt.title(r'Polariton Scattering Rates from $k_{0} \to k_{prime}$')
#--------------------------------------------------------------

temp = 50
k = np.linspace(0,1e6,1000)
k_prime = k
rate_mod = np.zeros((len(k),len(k_prime)), dtype = float)

for i in range(len(k)):
    for j in range(len(k_prime)):
        rate_mod[i,j] = psp.W_k_k_prime(k[i], k_prime[j], temp)
    print(i,'/',5000)
rate_mod = np.nan_to_num(rate_mod)
print('Fig1 Genned')
fig2 = plt.figure()
plt.pcolor(k*1e-6, k_prime*1e-6, rate_mod)
cbar2 = plt.colorbar()
cbar2.set_label(r'Scattering Rate ($s^{-1}$)')
plt.xlabel(r'|k| ($\mu m^{-1}$)')
plt.ylabel(r'$|k_{prime}|$ ($\mu m^{-1}$)')
plt.title(r'Polariton Scattering Rates from $|k| \to |k_{prime}|$')
print('Fig2 Genned')
plt.show()

