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
import emailer as email

hbar = 6.582119514e-16
kb = 8.6173303e-5
m0 = 9.10938356e-31
u = 4810
#--------------------------------------------------
num_sites = 50
num_steps = 1000
temp = 20
k_values = np.linspace(0,1000000,100)
gamma = 10e-6 / hbar #Polariton Decay Rate in the Cavity

class nLevelKMC:
    def __init__(self,num_sites,num_steps,temp,k_values):
        self.num_sites = num_sites
        self.num_steps = num_steps
        self.temp = temp
        self.k_values = k_values
        self.sites = np.zeros(self.num_sites, dtype = int)
        self.rate_lookup = np.zeros((len(self.k_values),len(self.k_values)), dtype = float)
        self.rates = np.zeros((len(self.k_values),len(self.k_values)), dtype = float)
        self.state_values = np.zeros(len(k_values), dtype = float)
        self.rate_sum = 0
        self.chosen_transition = 0
        self.time = 0
        self.time_array = [0]
        self.state_array = []
        self.total_momentum = []
        
    def initializeStates(self, state = 'gs'):
        if state == 'gs':
            for i in range(len(self.sites)):
                self.sites[i] = 0
        elif state == 'nonRes':
            for i in range(len(self.sites)):
                self.sites[i] = len(k_values) - 1
        elif state == 'random':
            for i in range(len(self.sites)):
                self.sites[i] = int(np.floor(len(k_values)*random()))
                
    def stateCounter(self):
        for i in range(len(self.k_values)):
            self.state_values[i] = (self.sites == i).sum()
                
    def makeRateLookup(self, decay = 0):
        for i in range(len(self.k_values)):
            for j in range(len(self.k_values)):
                self.rate_lookup[i,j]= psp.W_k_k_prime(self.k_values[i], self.k_values[j], self.temp)
    
    def updateRates(self, decay = 0):
        for i in range(len(self.k_values)):
            for j in range(len(self.k_values)):
                self.rates[i,j] = self.rate_lookup[i,j] * self.state_values[i] * (1 + self.state_values[j])
                if np.isnan(self.rates[i,j]):
                    self.rates[i,j] = 0
                
    def selectEvent(self):
        self.rate_sum = self.rates.sum()
        rand_rate = 0
        while rand_rate == 0:
            rand_rate = random()
        site = 0
        cumulative_rate = self.rates.flat[site]
        while cumulative_rate < (self.rate_sum * rand_rate):
            cumulative_rate += self.rates.flat[site + 1]
            site += 1
        self.chosen_transition = site
        self.chosen_transition = [int(site/len(k_values)), int(site%len(k_values))]
        
    def executeEvent(self):
        valid_sites = np.nonzero(np.in1d(self.sites,self.chosen_transition[0]))[0]
        chosen_site = choice(valid_sites)
        self.sites[chosen_site] = self.chosen_transition[1]
        
    def advanceTime(self):
        dt = 0
        rand_time = random()
        dt = -np.log(rand_time) / (self.rate_sum)
        self.time += dt
        self.time_array.append(self.time)
        
    def trackConfig(self):
        self.state_array.append(np.array(self.sites))
        self.total_momentum.append((max(k_values)/len(k_values)) * np.array(self.sites.sum()))
    
    def simulate(self):
        self.initializeStates()
        self.makeRateLookup()
        self.stateCounter()
        self.trackConfig()
        for i in range(self.num_steps):
            self.updateRates()
            self.selectEvent()
            self.executeEvent()
            self.advanceTime()
            if i < 10:
                self.nonResPump(1)#int(random()*10))
            self.stateCounter()
            self.trackConfig()

    def nonResPump(self, number):
        for i in range(number):
            self.sites = np.append(self.sites, (len(k_values))-1)
    
kmc = nLevelKMC(num_sites,num_steps,temp,k_values)
kmc.simulate()


    
e_loss = hbar*(u * kmc.total_momentum[-1])
ans = str('Total Phonon Energy Absorbed = '+ str(e_loss) +' eV')
print(ans)
fig1 = plt.figure()
plt.plot(kmc.time_array,kmc.total_momentum)
plt.title('Total Polariton Momentum of System wrt Time')
plt.xlabel('Time (Seconds)')
plt.ylabel(r'Total Momentum ($m^{-1}$)')

plt.show()

        














