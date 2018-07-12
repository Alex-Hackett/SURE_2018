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

hbar = 6.582119514e-16
kb = 8.6173303e-5
m0 = 9.10938356e-31
E_X_0 = 1.515
me = 0.067 * m0
mh = 0.45 * m0
a0 = 10e-9
n = 3.43
hbar_omega = 5e-3
De = -8.6
Dh = 5.7
u = 4810
rho = 5.3e3
lz = 1.192970175673805e-07
S = 1e-10
#V = lz * S
#V = 1e5
#V = 11546.4
V = 1000
c = 299792458

#Defining Dispersion Relation and Hopfield Coefficient
#Lower Polariton Only
def dis(k, extra = 0):
    E_C = ((hbar * c)/(n)) * (((np.pi)/(lz))**2 + (k**2))**(1/2)
    E_X = E_X_0 + 0*((hbar**2 * k**2)/(2*(me+mh)))
    
    E_L = 0.5 * (E_C + E_X - np.sqrt((E_C - E_X)**2 + (hbar_omega)**2))
    if extra:
        return E_L, E_C, E_X
    return E_L


#Exciton Hopfield Coefficient
def Hop_X(k):
    E_X = E_X_0 + ((hbar**2 * k**2)/(2*(me+mh)))
    E_L = dis(k)
    u_k = 1 / (np.sqrt(1 + ((hbar_omega/2)/(E_L - E_X))**2))
    return u_k
    
    

#Defining Helper Functions--------------------------------
def F(q):
    return (1 + ((q*a0)/(2))**2)**((-3)/(2))

def B(q):
    return ((8*np.pi**2)/(lz*q*((4*np.pi**2)-(lz**2 * q**2)))) * (np.sin((lz*q)/(2)))

def D(q):
    Fh = F((q*mh)/(me+mh))
    Fe = F((q*me)/(me+mh))
    return (De*Fh) - (Dh*Fe)

def d_kk(k,k_prime):
    E_k_prime = dis(k_prime)
    E_k = dis(k)
    delta_kk = (abs(E_k_prime - E_k))/(hbar*u)
    return delta_kk

def phonon_proj(k, k_prime):
    delta_kk = d_kk(k,k_prime)
    qz = ((delta_kk)**2 -(abs(k - k_prime))**2)**(0.5)
    return qz

def phonon_number(k, k_prime, temp):
    E_k_prime = dis(k_prime)
    E_k = dis(k)
    if E_k_prime - E_k > 0:
        return (np.exp((E_k_prime - E_k)/(kb*temp)) - 1)**(-1)
    elif E_k_prime - E_k < 0:
        return (np.exp((E_k_prime - E_k)/(kb*temp)) - 1)**(-1) + 1
    elif E_k_prime - E_k == 0:
        return 0
#Defining Main Scattering Rate Functions
def W_k_k_prime(k, k_prime, temp):
    '''
    E_X_0 = 1.515
    me = 0.067 * m0
    mh = 0.45 * m0
    a0 = 10e-9
    n = 3.43
    hbar_omega = 5e-3
    De = -8.6
    Dh = 5.7
    u = 4810
    rho = 5.3e3
    lz = 5e-9
    S = 1e-10
    V = lz * S
    V = 1e5
    V = 11546.4
    c = 299792458
    '''
    u_k = Hop_X(k)
    u_k_prime = Hop_X(k_prime)
    delta_k_k_prime = d_kk(k, k_prime)
    qz = phonon_proj(k, k_prime)
    first_term = (lz * (u_k*u_k_prime*delta_k_k_prime)**2)/(hbar * rho * V * qz * u**2)
    
    B_term = B(qz)
    B_sq = B_term**2
    D_term = D(abs(k - k_prime))
    D_sq = D_term**2
    
    N_ph = phonon_number(k, k_prime, temp)
    
    last_term = (delta_k_k_prime - abs(k - k_prime))
    
    W = first_term * B_sq * D_sq * N_ph * last_term
    return abs(W)


#Defining Parameters
'''
E_X_0 = 1.515
me = 0.067 * m0
mh = 0.45 * m0
a0 = 10e-9
n = 3.43
hbar_omega = 5e-3
De = -8.6
Dh = 5.7
u = 4810
rho = 5.3e3
lz = 5e-9
S = 1e-10
V = lz * S
V = 1e5
V = 11546.4
c = 299792458
#temp = 20
'''
def main():
    #Test
    '''
    k = np.linspace(0,10,1000)
    k_prime = np.linspace(0,10,1000)
    
    Ws = np.zeros((len(k), len(k_prime)),dtype = float)
    
    for i in range(len(k)):
        for j in range(len(k_prime)):
            W = W_k_k_prime(k[i],k_prime[j],temp)
            Ws[i,j] = W
        print(i,'/',len(k))
    
    fig1 = plt.figure()
    plt.pcolor(k, k_prime, Ws)
    plt.xlabel(r'k ($m^{-1}$)')
    plt.ylabel(r'$k_{prime}$ ($m^{-1}$)')
    plt.title('Polariton-Acoustic Phonon Scattering Rates')
    cbar = plt.colorbar()
    cbar.set_label('Scattering Rate ($s^{-1}$)')
    '''
    '''
    k = 0
    k_prime = np.linspace(-20,20,20000) * 1e4
    temp = 1
    Ws1 = np.zeros(len(k_prime), dtype = float)
    for i in range(len(k_prime)):
        Ws1[i] = W_k_k_prime(k, k_prime[i], temp)
        #print(i,'/',len(k_prime))
        
    
    k = 0
    k_prime = np.linspace(-20,20,20000) * 1e4
    temp = 5
    Ws5 = np.zeros(len(k_prime), dtype = float)
    for i in range(len(k_prime)):
        Ws5[i] = W_k_k_prime(k, k_prime[i], temp)
        #print(i,'/',len(k_prime))
        
        
        
    k = 0
    k_prime = np.linspace(-20,20,20000) * 1e4
    temp = 10
    Ws10 = np.zeros(len(k_prime), dtype = float)
    for i in range(len(k_prime)):
        Ws10[i] = W_k_k_prime(k, k_prime[i], temp)
        #print(i,'/',len(k_prime))
     '''   
        
        
    k = 0
    k_prime = np.linspace(-10e-6,10e6,20000)
    temp = 15
    Ws50 = np.zeros(len(k_prime), dtype = float)
    N_ph = np.zeros(len(k_prime), dtype = float)
    for i in range(len(k_prime)):
        Ws50[i], N_ph[i] = W_k_k_prime(k, k_prime[i], temp)
        
        #print(i,'/',len(k_prime))
    
    
    fig2 = plt.figure()
    #plt.plot(k_prime, Ws1, label = '1 Kelvin')
    #plt.plot(k_prime, Ws5, label = '5 Kelvin')
    #plt.plot(k_prime, Ws10, label = '10 Kelvin')
    plt.plot(k_prime, Ws50, label = '15 Kelvin')
    #plt.plot(k_prime, N_ph, label = 'Phonon Number')
    plt.xlabel(r'$k_{prime}$ $m^{-1}$')
    plt.ylabel(r'Rate ($s^{-1}$)')
    plt.title(r'Scattering Rate $k_0 \to k_{prime}$')
    plt.legend()
    
    
    
    fig3 = plt.figure()
    k = np.linspace(-10e6,10e6,10000) 
    lp, ec, ex = dis(k,1)
    plt.plot(k,lp, label = 'Lower Polariton Branch')
    plt.plot(k,ec, label = 'Cavity Photon Dispersion')
    plt.plot(k,ex, label = 'Exciton Dispersion')
    plt.xlabel(r'k ($m^{-1}$)')
    plt.ylabel(r'Energy (eV)')
    plt.legend()
    
    '''
    temp = 20
    k = np.linspace(-20,20,100) * 1e4
    k_prime = np.linspace(-20,20,100) * 1e4
    W = np.zeros((len(k), len(k_prime)), dtype = float)
    for i in range(len(k)):
        for j in range(len(k_prime)):
            W[i,j] = W_k_k_prime(k[i], k[j],temp)
    fig5 = plt.figure()
    plt.pcolor(k, k_prime, W)
    plt.xlabel(r'k ($m^{-1}$)')
    plt.ylabel(r'$k_{prime}$ ($m^{-1}$)')
    plt.title('Polariton-Acoustic Phonon Scattering Rates')
    cbar = plt.colorbar()
    cbar.set_label('Scattering Rate ($s^{-1}$)')
    '''   
    
    
    
    plt.show()
if __name__ == '__main__':
    main()
    
    
