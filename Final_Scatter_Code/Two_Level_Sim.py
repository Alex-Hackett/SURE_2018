# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 12:51:31 2018

@author: Alex
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
import polariton_scattering_package as psp

hbar = 6.582119514e-16
kb = 8.6173303e-5
m0 = 9.10938356e-31
#--------------------------------------------------
num_sites = 1000
num_steps = 6520000
temp = 20
k0 = 0
k1 = 5000



class kineticMonteCarlo:
    def __init__(self,num_sites,num_steps,temp):
        self.num_sites = num_sites
        self.num_steps = num_steps
        self.sites = np.zeros(self.num_sites)
        self.rates = np.zeros(self.num_sites)
        self.k1s = []
        self.k0s = []
        self.k1count = 0
        self.k0count = 0
        self.time = 0
        self.selected_site = 0
        self.rate_sum = 0
        self.k1_to_k0 = psp.W_k_k_prime(k1,k0,temp)
        self.k0_to_k1 = psp.W_k_k_prime(k0,k1,temp)
        self.time_array = []
        self.hack_count = 0
    def make_initial_states(self, gs = 0):
        if gs:
            for i in range(self.num_sites):
                self.sites[i] = 1
        elif not gs:
            for i in range(self.num_sites):
                self.sites[i] = choice([0,1])
        self.sites[4] = 1
        
    def count_states(self):
        self.k1count = (self.sites == 1).sum()
        self.k0count = (self.sites == 0).sum()
        self.k1s.append(self.k1count)
        self.k0s.append(self.k0count)
        
    def sum_rates(self):
        self.rate_sum = 0
        for i in range(self.num_sites):
            self.rate_sum += self.rates[i]
        
    def choose_event(self):
        rand_rate = 0
        while rand_rate == 0:
            rand_rate = random()
        x = 0
        rate_cum = self.rates[x]
        while rate_cum < self.rate_sum * rand_rate:
            rate_cum += self.rates[x+1]
            x += 1
        self.selected_site = x
        
    def update_rates(self):
        k1_k0 = self.k1_to_k0 * (self.k1count) * (1+self.k1count)
        k0_k1 = self.k0_to_k1 * (self.k0count) * (self.k1count)
        for i in range(self.num_sites):
            if self.sites[i] == 1:
                self.rates[i] = abs(k1_k0)
            elif self.sites[i] == 0:
                self.rates[i] = abs(k0_k1)
    
    def do_event(self):
        if self.sites[self.selected_site] == 0:
            self.sites[self.selected_site] = 1
        elif self.sites[self.selected_site] == 1:
            self.sites[self.selected_site] = 0
            
    def advance_time(self):
        if self.rate_sum == 1:
            print('Two+ Problems in a Row!')
        if self.rate_sum == 0:
            self.hack_count += 1
            self.rate_sum =1
            print(self.hack_count,'th Divide by Zero Hack in Time Advancment')
        dt = 0
        rand_time = random()
        dt = -np.log(rand_time) / (self.rate_sum)
        self.time += dt
        self.time_array.append(self.time)
        
    def simulate(self):
        self.make_initial_states()
        for i in range(self.num_steps):
            self.count_states()
            self.update_rates()
            self.sum_rates()
            self.choose_event()
            self.do_event()
            self.advance_time()
            #print(i,'/',self.num_steps)
            

kmc = kineticMonteCarlo(num_sites, num_steps, temp)
kmc.simulate()
fig1 = plt.figure()
plt.plot(kmc.time_array,kmc.k1s, label = 'Number of k1 States')
plt.plot(kmc.time_array,kmc.k0s, label = 'Number of k0 States')
plt.xlabel('Time (Seconds)')
plt.ylabel('Number of States')
plt.title('Binary Kinetic Monte Carlo at 20K')
plt.legend()
plt.show()


















