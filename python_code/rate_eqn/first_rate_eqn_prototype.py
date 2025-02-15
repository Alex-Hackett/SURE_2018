#!/bin/python

'''
Just a quick attempt to describe the acoustic phonon-polariton scattering rates
'''

#Imports
#Imports
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import astropy as ap
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii


#Defining Variables
#All masses in free electron mass
#All energies in eV
#All distances in nm, momenta in 1/nm

num_pol = 1e6
L = 10.#Nanowire Size
R = 1.#Nanowire Radius
V = np.pi * R*R * L 
rho = 5318 * 1e-67#GaAs density
u = 3350 * 1e9 #Speed of sounds
hbar = 6.582119514e-16#Planck's constant
omega_r = 10e-3#Rabi Splitting
boltzmann = 8.6173303e-5 #Boltzmann's constant
a_e = -7 #Phonon lattice deformation at electrons
a_h = 2.7 #Phonon lattice deformation at holes
a_B = 10 #Exciton Bohr radius
m_e = 0.067 #Effective electron mass
m_h = 0.18 #Effective hole mass
epsilon = 11. #Dielectric constant for GaAs


#Simulation parameters
k1 = 1/0.01#Initial polariton wavevector
k2 = 1/0.02#Final polariton wavevector
dk = k1 - k2
qz = 0.005#Phonon vertical wavevector
beta = boltzmann * 273


#Polariton Dispersion Relation
E_p1 = np.sqrt(u**2 * k1**2 * 1./epsilon)
E_p2 = np.sqrt(u**2 * k2**2 * 1./epsilon)

#Excitonic Fraction function
X_p1 = 2./ np.sqrt(4 + (abs(E_p1)/omega_r)**2)
X_p2 = 2./ np.sqrt(4 + (abs(E_p2)/omega_r)**2)

#Overlap integrals
I_par_e = (1 + (((m_e)/(m_e + m_h)) * abs(dk)*a_B)**2)**(-3/2)
I_par_h = (1 + (((m_h)/(m_e + m_h)) * abs(dk)*a_B)**2)**(-3/2)

I_perp_e = ((np.pi)**2) * (np.sin((qz * L)/2.)) * 1/((0.5 * qz *L) * ((np.pi)**2 - (0.5 * qz * L)**2))

I_perp_h = ((np.pi)**2) * (np.sin((qz * L)/2.)) * 1/((0.5 * qz *L) * ((np.pi)**2 - (0.5 * qz * L)**2))

#Phonon Assisted Scattering Rate
W = (L/(rho * u * V)) * ((abs(dk)**2 + qz**2)/(abs(hbar * u * qz))) * (abs(X_p1) * abs(X_p2))**2 * ((a_e*I_par_e*I_perp_e) - (a_h*I_par_h*I_perp_h))**2


#Test rate from 1/0.2 to 1/0.1
P_dot = W * num_pol * (1 + 1/(np.exp(hbar * qz * u * beta) + 1))



t = 0
t_ar = []
t_ar.append(t)
p_ar = []
p_ar.append(num_pol)
while t < 1000:
    print('Polariton Number = ',num_pol)
    num_pol = num_pol - (W * num_pol * (1 + 1/(np.exp(hbar * qz * u * beta) + 1)))
    p_ar.append(num_pol)
    t += 0.01
    t_ar.append(t)

fig1 = plt.plot(t_ar, p_ar)
plt.show()



