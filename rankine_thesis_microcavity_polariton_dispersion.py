# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:29:04 2018

@author: Alex

Function module for getting the lower and upper branch polarition
dispersions as a function of energy from the in plane wavevector 
and other parameters
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
m_h = 0.18 * 9.10938356e-31 #Hole effective mass
m_e = 0.067 * 9.10938356e-31 #Electron effective mass
hbar = 6.582119514e-16 #in eV.s
lz = 300e-9 #Cavity Lenght in Meters
M = 1 #Mode of coupled Light
c = 299792458 
n0 = 3.9476 #Refective index of CAVITY
qz = 2 * np.pi * M * (1/lz)
omega_ex_0 = 3.54
omega_cav_0 = omega_ex_0
#omega_ex_0 = omega_ex_0/hbar
rabi = 50e-3#Rabi splitting in eV
#rabi = rabi/hbar
little_m = 0.5e-4 * 9.10938356e-31
big_m = 0.25 * 9.10938356e-31
def PDis(mod_k, lz, M, n0, rabi, omega_ex_0,upflag = 0):
    '''
    Parameters are 
    mod_k, the magnitude of the in plane wavevector
    upflag, wip atm
    lz, the effective optical resonator length
    M, the photon mode (probably unity)
    n0, the refractive index of the resonator/cavity
    rabi, the Rabi splitting of the material in eV
    omega_ex_0, the zero wavevector exciton frequency, in eV
    
    Outputs are omega_l*hbar and omega_u*hbar, the lower and upper branch
    polariton energies in eV
    '''
    hbar = 6.582119514e-16 #in eV.s
    
    c = 299792458 
    #n0 = 3.9476 #Refective index of CAVITY
    #qz = 2 * np.pi * M * (1/lz)
    #k_perp = n0 * ((2*np.pi)/(lz))
    
    #omega_cav = hbar*c*(1/n0) * np.sqrt(mod_k**2 + k_perp**2)
    omega_cav = omega_cav_0 + ((hbar**2 * mod_k**2)/(2*little_m))
    
    omega_ex = omega_ex_0 + ((hbar**2 * mod_k**2)/(2*(big_m)))
    

    #Lower Branch Polariton Dispersion
    omega_l = 0.5 * (omega_ex+omega_cav - np.sqrt(rabi**2 + (omega_ex-omega_cav)**2))
    #Upper Branch Polariton Dispersion 
    omega_u = 0.5 * (omega_ex+omega_cav + np.sqrt(rabi**2 + (omega_ex-omega_cav)**2))
    return omega_l, omega_u, omega_cav, omega_ex

def main():
    k_array = np.linspace(-2,2,100) * 1e-6
    omega_array_l = np.zeros(len(k_array))
    omega_array_u = np.zeros(len(k_array))
    cavity_array = np.zeros(len(k_array))
    exciton_array = np.zeros(len(k_array))
    
    j = 0
    for i in k_array:
        omega_array_l[j], omega_array_u[j],cavity_array[j],exciton_array[j] = PDis(i,lz,1,n0,rabi,omega_ex_0)
        
        j += 1
    fig1 = plt.figure()
    plt.plot(k_array*1e6,omega_array_l , label = 'Lower Polariton Branch')
    plt.plot(k_array*1e6,omega_array_u , label = 'Upper Polarition Branch')
    plt.plot(k_array*1e6,cavity_array , label = 'Photon Dispersion')
    plt.plot(k_array*1e6,exciton_array , label = 'Exciton Dispersion')
    plt.xlabel(r'Wavevector |K| ($\mu m^{-1}$)')
    plt.ylabel(r'E (eV)')
    #plt.ylim(1600,1700)
    plt.legend()
    plt.title('Polariton Dispersion Curves')
    plt.show()
    
if __name__ == '__main__':
    main()


    



