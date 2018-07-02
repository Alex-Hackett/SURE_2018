# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:29:04 2018

@author: Alex

Function module for getting the lower and upper branch polarition
dispersions as a function of energy from the in plane wavevector 
and other parameters
and the excitionic fraction
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


def PDis(mod_k, upflag = 0, lz, M, n0, rabi, omega_ex_0):
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
    #lz = 300e-6 #Cavity Lenght in Meters
    #M = 1 #Mode of coupled Light
    c = 299792458 
    #n0 = 3.9476 #Refective index of CAVITY
    qz = 2 * np.pi * M * (1/lz)
    #omega_ex_0 = 1.557/hbar
    omega_ex_0 = omega_ex_0/hbar
    #rabi = 30 * 1e-3#Rabi splitting in eV
    rabi = rabi/hbar
    
    omega_cav = c/n0 * np.sqrt(qz**2 + mod_k**2) #Cavity Mode dispersion
    omega_ex = omega_ex_0 + (1e-4 * mod_k**2) #Exciton Dispersion
    
    #Lower Branch Polariton Dispersion
    omega_l = 0.5 * ((omega_cav + omega_ex) - np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
    #Upper Branch Polariton Dispersion 
    omega_u = 0.5 * ((omega_cav + omega_ex) + np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
    return omega_l*hbar, omega_u*hbar

def main():
    k_array = np.linspace(-10,10,10000) * 1e6
    omega_array_l = np.zeros(len(k_array))
    omega_array_u = np.zeros(len(k_array))
    j = 0
    for i in k_array:
        omega_array_l[j], omega_array_u[j] = PDis(i)
        
        j += 1
    fig1 = plt.figure()
    plt.plot(k_array*1e-6,omega_array_l * 1e3, label = 'Lower Polariton Branch')
    #plt.plot(k_array*1e-9,omega_array_u * 1e3, label = 'Upper Polarition Branch')
    plt.xlabel(r'Wavevector |K| ($\mu m^{-1}$)')
    plt.ylabel(r'E (meV)')
    plt.legend()
    plt.title('Polariton Dispersion Curves')
    plt.show()
    
if __name__ == '__main__':
    main()


    



