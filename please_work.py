# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:29:04 2018

@author: Alex

Function module for getting the lower and upper branch polarition
dispersions as a function of energy from the in plane wavevector 
and other parameters
and the excitionic fraction
"""
import sympy as sy
import scipy as sp
import os
import astropy as ap
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import math
import numpy as np
import scipy.integrate as integrate
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import sys
import time
from random import randrange, random, choice
import os.path
omega_ex_0 = 1417.5e-3
omega_cav_0 = omega_ex_0
#omega_ex_0 = omega_ex_0/hbar
rabi = 3.5 * 1e-3#Rabi splitting in eV
#rabi = rabi/hbar
little_m = 2.5e-5 * 9.10938356e-31
big_m = 0.25 * 9.10938356e-31
hbar = 6.582119514e-16
import rankine_thesis_microcavity_polariton_dispersion as newdis


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
    return E_l#, E_p, e_cav, e_ex















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


def scatter_rate(qz,omega_ex_0, rabi, n0, k1, k2, L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V):
    '''
    This is the main scattering function
    TODO, Monday, finish this
    TODO, Tuesday. test this, attempt to obtain results from kinetic 
    MC paper
    '''
    
    
    d_k = abs(k1 - k2)
    first_term = L**2 / (rho * u * V)
    second_term = ((abs(d_k))**2 + qz**2)/(abs(hbar * u * qz))
    pol_E_k1 = polariton_dispersion(k1)
    pol_E_k2 = polariton_dispersion(k2)
    X_1 = XFrac(rabi, pol_E_k1)
    X_2 = XFrac(rabi, pol_E_k2)
    exciton_term = abs(X_1 * X_2)**2
    perp_e = IPerp(qz, lz)
    par_e = IPar(d_k, m_e, m_h, a_b, 'e')
    perp_h = perp_e
    par_h = par_e = IPar(d_k, m_e, m_h, a_b, 'h')
    
    integral_term = ((a_e * par_e * perp_e) - (a_h * par_h * perp_h))**2
    
    return first_term * second_term * exciton_term * integral_term

def makeqz(theta,k1,k2):
    d_k_par_sq = k1**2 +k2**2 - (2*k1*k2*np.cos(theta))
    qz = np.sqrt((k1-k2)**2 - d_k_par_sq)
    return qz

def maketheta_max(k1,k2):
    cos_t_max  = (k1**2 + k2**2 - (k1-k2)**2)/(2*k1*k2)
    if cos_t_max > 1:
        cos_t_max = 1
    if cos_t_max < -1:
        cos_t_max = -1
    return (cos_t_max)
    



#Defining Constants according to KMC paper

omega_ex_0 = (1.557) #Zero Momentum Exciton Freq (eV)
rabi = (10 * 1e-3) #Rabi Splitting (eV)
n0 = 3.857 #Refractive Index of GaAs


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
k2 = np.linspace(-10,10,500) * 1e-9 #Incoming Wavevector
#k2 = np.zeros(len(k2))
k1 = np.linspace(-10,10,500) * 1e-9 #Outgoing wavevector
rates = np.zeros((len(k1),len(k2)), dtype = float) #Array to hold rates


for i in range(len(k2)):
    for j in range(len(k1)):
        #qz = abs(abs(k1[j]) - abs(k2[i]))
        
        #max_theta = maketheta_max(k1[j],k2[i])
        omega_k2 = polariton_dispersion(k2[i])
        omega_k1 = polariton_dispersion(k1[j])
        qz = np.sqrt(abs(((abs(omega_k1) - abs(omega_k2))/(hbar*u))**2 - (abs(k1[j]) - abs(k2[i]))**2))
        rate = scatter_rate(qz, omega_ex_0, rabi, n0, abs(k1[j]), abs(k2[i]), L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V)
        #rate, error = integrate.quad(scatter_rate, -abs(abs(k1[j])-abs(k2[i])), abs(abs(k1[j])-abs(k2[i])), args = (omega_ex_0, rabi, n0, (k1[j]), (k2[i]), L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V), limit = 2000)
        if np.isnan(rate):
            rate = 0
        else:
            rates[j,i] = rate * hbar
        #print(i,'/',len(k2) - 1)
        #print(j,'/',len(k1))
    print(i,'/',len(k2) - 1)



fig1 = plt.figure()
plt.pcolor(k2 , k1, ((rates * 1e9)))#/(max(max(x) for x in rates)) * 256))
plt.xlabel(r'$k_{2}$ ($\mu m^{-1}$)')
plt.ylabel(r'$k_{1}$ ($\mu m^{-1}$)')
#plt.xlim(-10e6,10e6)
#plt.ylim(-10e6,10e6)
plt.title('Phonon-Polariton Scattering Rates')
cbar = plt.colorbar()
cbar.set_label('Scattering Rate (neV)')
plt.show()
'''
fig1 = plt.figure()
plt.plot(k1*1e-6, rates[:,0]*1e9)
plt.xlabel(r'$k_{1}$ ($\mu m^{-1}$)')
plt.ylabel(r'Scattering Rate (neV)')
plt.title('Polariton-Phonon Scattering Rate to $|k|=0$ as a function of Wavevector')
'''
'''
plt.show()
    
    
E1,dum = PDis(k1[j], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        E2,dum = PDis(k2[i], lz, 1, n0, rabi, omega_ex_0, upflag = 0)
        qz = (((E2-E1)/(hbar*u))**2 - (k1[j]-k2[i])**2)**0.5

'''


