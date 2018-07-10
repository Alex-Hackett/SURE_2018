# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:29:04 2018

@author: Alex

Function module for getting the lower and upper branch polarition
dispersions as a function of energy from the in plane wavevector 
and other parameters
and the excitionic fraction
"""
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

hbar = 6.582119514e-16
kb = 8.6173303e-5







def PDis(mod_k, lz, M, n0, rabi, omega_ex_0, upflag = 0):
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
    
    omega_cav = (c/n0) * np.sqrt(qz**2 + mod_k**2) #Cavity Mode dispersion
    omega_ex = omega_ex_0 + (1e-4 * mod_k**2) #Exciton Dispersion
    
    #Lower Branch Polariton Dispersion
    omega_l = 0.5 * ((omega_cav + omega_ex) - np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
    #Upper Branch Polariton Dispersion 
    omega_u = 0.5 * ((omega_cav + omega_ex) + np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
    return omega_l*hbar, omega_u*hbar







'''
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
'''



def XFrac(rabi, pol_E):
    '''
    Simple function returns the excitionic fraction as a function of 
    the rabi splitting and the polariton energy at some K (use the polarition
    dispersion function to acquire this)
    '''
    
    return 2/(np.sqrt(4 + (((pol_E)/(rabi)))**2))


def IPar(delta_k, m_e, m_h, a_b, e_or_h):
    '''
    This is the parallel overlap intergral
    '''
    if str(e_or_h) == 'e':
        return (1 + (((m_e)/(m_h + m_e))*abs(delta_k)*(a_b))**2)**(-3/2)
    elif str(e_or_h) == 'h':
        return (1 + (((m_h)/(m_h + m_e))*abs(delta_k)*(a_b))**2)**(-3/2)
    

def IPerp(qz, lz):
    '''
    gets the perpiducular overlap integral with the phonon wavefunction
    should be independent of whether or not the interacting object is an
    electron or hole
    '''
    return (((np.pi)**2) * np.sin((qz*lz)/(2)))/(((qz*lz)/(2)) * (((np.pi)**2) - ((qz*lz)/(2))**2))


def scatter_rate(omega_ex_0, rabi, n0, k1, k2, L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V):
    '''
    This is the main scattering function
    TODO, Monday, finish this
    TODO, Tuesday. test this, attempt to obtain results from kinetic 
    MC paper
    '''
    def PDisdiffy(mod_k):
        omega_ex_0 = (1.557) #Zero Momentum Exciton Freq (eV)
        rabi = (10 * 1e-3) #Rabi Splitting (eV)

        #mod_k = mod_k**2
        hbar = 6.582119514e-16 #in eV.s
        #lz = 300e-6 #Cavity Lenght in Meters
        #M = 1 #Mode of coupled Light
        c = 299792458 
        #n0 = 3.9476 #Refective index of CAVITY
        qz = 2 * np.pi * 1 * (1/lz)
        #omega_ex_0 = 1.557/hbar
        omega_ex_0 = omega_ex_0/hbar
        #rabi = 30 * 1e-3#Rabi splitting in eV
        rabi = rabi/hbar
        
        omega_cav = (c/n0) * np.sqrt(qz**2 + mod_k**2) #Cavity Mode dispersion
        omega_ex = omega_ex_0 + (1e-4 * mod_k**2) #Exciton Dispersion
        
        #Lower Branch Polariton Dispersion
        omega_l = 0.5 * ((omega_cav + omega_ex) - np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
        #Upper Branch Polariton Dispersion 
        omega_u = 0.5 * ((omega_cav + omega_ex) + np.sqrt((omega_cav - omega_ex)**2 + (4*rabi**2)))
        if k2 == 0:
            return np.inf
        return omega_l*hbar #, omega_u*hbar
    
    
    
    
    d_k = abs(abs(k1) - abs(k2))
    d_e = (k1-k2) * hbar * u
    first_term = ((2*np.pi)/(hbar)) * (S/(2*np.pi)**2)
    diff_2_E = spmisc.derivative(PDisdiffy,(k2), dx = 1e-9)  * (1/(2*k2))
    diff_2_E = 1
    second_term = (d_e)/(2*diff_2_E)
    pol_E_k1, dum = PDis(k1, lz, 1, n0, rabi, omega_ex_0)
    pol_E_k2, dum = PDis(k2, lz, 1, n0, rabi, omega_ex_0)
    X_1 = XFrac(rabi, pol_E_k1)
    X_2 = XFrac(rabi, pol_E_k2)
    
    
    exciton_term = abs(X_1 * X_2)
    temp = 300
    
    third_term = ((hbar*d_k)/(2*rho*S*u))
    N_ph = 0
    if d_e > 0:
        N_ph = (np.exp(d_e/(kb*temp)) - 1)**-1
    elif d_e < 0:
        N_ph = ((np.exp(d_e/(kb*temp)) - 1)**-1) + 1
        
    
    
    
    par_e = IPar(d_k, m_e, m_h, a_b, 'e')
   
    par_h = par_e = IPar(d_k, m_e, m_h, a_b, 'h')
    
    integral_term = ((a_e * par_e)- (a_h * par_h))**2
    
    #qz = np.sqrt(d_k**2 - (k1**2 + k2**2 - (2*k1*k2*np.cos(theta))))
    if (2*k1*k2) == 0:
        c = np.inf
    else:
        c = (k1**2 + k2**2 - d_k**2)/(2*k1*k2)
    if c>1:
        cos_theta_max = 1
    elif c<-1:
        cos_theta_max = -1
    else:
        cos_theta_max = c
        
    theta_max = np.arccos(cos_theta_max)
    
    def r_integral(theta):
       return ((2*d_k)/(hbar * u * (k1**2 + k2**2 - (2*k1*k2*np.cos(theta))))) * ((np.sin((lz*(k1**2 + k2**2 - (2*k1*k2*np.cos(theta))))/2) * np.pi**2)/((0.5*lz*(k1**2 + k2**2 - (2*k1*k2*np.cos(theta))))*(np.pi**2 - (0.5*lz*(k1**2 + k2**2 - (2*k1*k2*np.cos(theta))))**2)))**2
    
    R_prime, err = integrate.quad(r_integral, 0, cos_theta_max, limit = 5000)
    R_prime = R_prime * 2 * integral_term
    
    return abs(first_term * second_term * exciton_term * third_term * R_prime * N_ph)


#Defining Constants according to KMC paper
    
omega_ex_0 = (1.557) #Zero Momentum Exciton Freq (eV)
rabi = (10 * 1e-3) #Rabi Splitting (eV)
n0 = 3.857 #Refractive Index of GaAs
S = 8e-6 * 10e-9#Quantization Surface

L = 8e-6 #Microcavity Length
lz = 10e-9 #Quantum Well Width
u = 3350 #Speed of sound in the cavity
rho = 5318 #Density of GaAs
m_h = 0.18 * 9.10938356e-31 #Hole effective mass
m_e = 0.067 * 9.10938356e-31 #Electron effective mass
a_b = 10e-9 #Exciton Bohr radius
a_h = 2.7 #Hole Lattice deformation constant (eV)
a_e = -7 #Electron lattice deformation constant (eV)
V = np.pi * lz * (L)**2 #QW Effective volume
k2 = np.linspace(-10,10,101) * 1e6 #Incoming Wavevector
#k2 = np.zeros(len(k2))
k1 = np.linspace(-10,10,101) * 1e6 #Outgoing wavevector
rates = np.zeros((len(k1),len(k2)), dtype = float) #Array to hold rates


for i in range(len(k2)):
    for j in range(len(k1)):
        #omega_k2,dum = PDis(k2[i], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        #omega_k1,dum = PDis(k1[j], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        #qz = np.sqrt(abs(((abs(omega_k1) - abs(omega_k2))/(u))**2 - (abs(k1[j]) - abs(k2[i]))**2))
        rate = scatter_rate(omega_ex_0, rabi, n0, abs(k1[j]), abs(k2[i]), L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V)
        #rate, error = integrate.quad(scatter_rate, -abs(k1[j] - k2[i]), abs(k1[j] - k2[i]), args = (omega_ex_0, rabi, n0, abs(k1[j]), abs(k2[i]), L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V))
        if np.isnan(rate):
            rate = 0
        else:
            rates[j,i] = rate * hbar
        #print(i,'/',len(k2) - 1)
        #print(j,'/',len(k1))
    print(i,'/',len(k2) - 1)

fig1 = plt.figure()
plt.pcolor(k2, k1, ((rates)))#/(max(max(x) for x in rates)) * 256))
plt.xlabel(r'$k_{2}$ ($m^{-1}$)')
plt.ylabel(r'$k_{1}$ ($m^{-1}$)')
#plt.xlim(-10e6,10e6)
#plt.ylim(-10e6,10e6)
plt.title('Phonon-Polariton Scattering Rates')
cbar = plt.colorbar()
cbar.set_label('Scattering Rate (eV)')
'''
fig1 = plt.figure()
plt.plot(k1*1e-6, rates[:,0]*1e9)
plt.xlabel(r'$k_{1}$ ($\mu m^{-1}$)')
plt.ylabel(r'Scattering Rate (neV)')
plt.title('Polariton-Phonon Scattering Rate to $|k|=0$ as a function of Wavevector')
'''

plt.show()
    
    
'''
E1,dum = PDis(k1[j], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        E2,dum = PDis(k2[i], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        qz = (((E2-E1)/(hbar*u))**2 - (k1[j]-k2[i])**2)**0.5
'''





