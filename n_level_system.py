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
#--------------------------------------------------
num_sites = 100
num_steps = 100
temp = 20
k_values = np.linspace(0,100000,10)

class nLevelKMC:
    def __init__(self,num_sites,num_steps,temp,k_values):
        self.num_sites = num_sites
        self.num_steps = num_steps
        self.temp = temp
        self.k_values = k_values
        self.sites = np.zeros(self.num_sites)
        self.rate_lookup = np.zeros((len(self.k_values),len(self.k_values)), dtype = float)
        self.rates = np.zeros((len(self.k_values),len(self.k_values)), dtype = float)
        self.state_values = np.zeros(len(k_values), dtype = float)
        
    def initializeStates(self, state = 'gs'):
        if state == 'gs':
            for i in range(self.sites):
                self.sites[i] = 0
        elif state == 'nonRes':
            for i in range(self.sites):
                self.sites[i] = 99
        elif state == 'random':
            for i in range(self.sites):
                self.sites[i] = int(np.floor(100*random()))
                
    def stateCounter(self):
        for i in range(len(self.k_values)):
            self.state_values[i] = (self.sites == i).sum()
                
    def makeRateLookup(self, decay = 0):
        for i in len(self.k_values):
            for j in len(self.k_values):
                self.rate_lookup[i,j]= psp.W_k_k_prime(self.k_values[i], self.k_values[j], self.temp)
    
    def updateRates(self, decay = 0):
        for i in len(self.k_values):
            for j in len(self.k_values):

            
        














