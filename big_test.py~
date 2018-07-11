# -*- coding: utf-8 -*-
"""
Test code for binary state kinetic monte carlo utilizing realistic scattering
rates from polarition_scatter_pack.py
"""
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

import polarition_scatter_pack as psp
kb = 8.6173303e-5

#Constants
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

hbar = 6.582119514e-16
#Select state k as 10um, state 0 as 0um
k1 = 10e6
k2 = 0
class Kinetic_Monte_Carlo:
    def __init__(self, number_sites, number_steps, temp):
        self.num_sites = number_sites
        self.num_steps = number_steps
        self.temp = temp
        self.rate_k_to_0 = 0
        self.rate_0_to_k = 0
        self.site = np.zeros(self.num_sites, dtype = float)
        self.mcs = 0
        self.rate = np.zeros(self.num_sites, dtype = float)
        self.rate_sum = 0
        self.time = 0
        self.rand_rate = 0.
        self.rand_time = 0.
        self.selected_site = 0
        self.delta_t = 0
        self.N_ph_k = ((np.exp(((k1*u*hbar)/(self.temp * kb))))**-1)
        self.k_num = 0
        self.zero_num = 0
        self.k_array = []
        self.zero_array = []
        self.time_array = []
    def initialize_system(self):
        for i in range(self.num_sites):
            self.site[i] = 0
            self.update_rates()
            
            
            
    def update_rates(self):
        
        qz = k1
        k_count = (self.site > 0).sum()
        zero_count = (self.site == 0).sum()
        
        #self.N_ph_k = ((np.exp(((k1*u*hbar)/(self.temp * kb))))**-1)
        
        W_zero_k = psp.scatter_rate(qz, omega_ex_0, rabi, n0, k2, k1, L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V)
        W_zero_k = W_zero_k * (self.N_ph_k)
        
        
        W_k_zero = psp.scatter_rate(qz, omega_ex_0, rabi, n0, k1, k2, L, lz, u, rho, m_e, m_h, a_b, a_e, a_h, V)
        W_k_zero = W_k_zero * (1 + self.N_ph_k)
        
        
        self.rate_0_to_k = (-k_count * W_k_zero * (1 + zero_count)) + ((1 + k_count) * W_zero_k * zero_count)
        self.rate_k_to_0 = (k_count * W_k_zero * (1 + zero_count)) - ((1 + k_count) * W_zero_k * zero_count)
        if self.rate_0_to_k < 0:
            self.rate_0_to_k = 0
            
        if self.rate_k_to_0 < 0:
            self.rate_k_to_0 = 0
        
        for i in range(self.num_sites):
            if self.site[i] == 0:
                    self.rate[i] = self.rate_0_to_k
            elif self.site[i] > 0:
                    self.rate[i] = self.rate_k_to_0
                 
    def count_states(self):
        self.k_num = 0
        self.zero_num = 0
        self.k_num = (self.site > 0).sum()
        self.zero_num = (self.site == 0).sum()
        self.k_array.append(self.k_num)
        self.zero_array.append(self.zero_num)

    def compute_rate_sum(self):
        self.rate_sum = 0
        for i in range(self.num_sites):
            self.rate_sum += self.rate[i]
            
    def choose_site(self):
        self.rand_rate = 0
        while self.rand_rate == 0:
            self.rand_rate = random()
        x = 0
        rate_cum = self.rate[x]
        while rate_cum < self.rate_sum * self.rand_rate:
            rate_cum += self.rate[x+1]
            x += 1
        self.selected_site = x
        
    def do_event(self):
        if self.site[self.selected_site] == 0:
            self.site[self.selected_site] = k1
            self.N_ph_k = self.N_ph_k - 1
        elif self.site[self.selected_site] > 0:
            self.site[self.selected_site] = 0
            self.N_ph_k += 1
        self.update_rates()
        
    def advance_time(self):
        self.dt = 0
        self.rand_time = random()
        self.dt = -np.log(self.rand_time) / (self.rate_sum)
        self.time += self.dt
        self.time_array.append(self.time)
    def evolve(self):
        for self.mcs in range(self.num_steps):
            self.compute_rate_sum()
            self.choose_site()
            self.do_event()
            self.advance_time()
            self.update_rates()
            self.count_states()
            self.mcs += 1
            print('Monte Carlo Step ',self.mcs,' of ',self.num_steps)
            
    def simulate(self):
        self.initialize_system()
        self.evolve()
        data = Table([self.time_array, self.k_array,self.zero_array], names = ['Time','k_count','zero_count'])
        filename = 'big_old_allnight_test.txt'
        ascii.write(data,filename)
        self.vis()
        
    def vis(self):
        fig1 = plt.figure()
        plt.plot(self.time_array, self.k_array, label = r'Number of k = 10$\mu m^{-1}$ states')
        plt.plot(self.time_array, self.zero_array, label = r'Number of k = 0$\mu m^{-1}$ states')
        plt.legend()
        plt.xlabel('Time (Seconds)')
        plt.ylabel('Number of States')
        plt.show()
        
        

kmc = Kinetic_Monte_Carlo(100,int(1e8),300)
kmc.simulate()

    

        
            
    
        
        
        
            
        


