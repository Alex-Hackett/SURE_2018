# -*- coding: utf-8 -*-
"""
This is a translation of the example KMC C code found in 
doi:10.1016/j.cma.2008.03.010 into python 3.6
A is state K, B is state 0
"""

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

NUM_SITES = 64
NUM_STEPS = 64000






class KMC:
    def __init__(self, sites_num, steps_num):
        self.RATE_A2B = 0.
        self.RATE_B2A = 0.
        self.temp = 100
        self.num_sites = sites_num
        self.num_steps = steps_num
        self.site = np.zeros(sites_num, dtype = int)
        self.prod = np.zeros(sites_num, dtype = int)
        self.mcs = 0
        self.rate = np.zeros(sites_num, dtype = float)
        self.rate_sum = 0.
        self.t = 0.
        self.pct_B_avg = 0.
        self.rand_rate = 0.
        self.rand_time = 0.
        self.selected_site = 0
        self.pct_B = 0.
        self.dt = 0.
        self.rate_listA2B = []
        self.rate_listB2A = []
        self.B_list = [0]
        self.g_pulse = True
        self.gamma_k =  1e12 * 2.6364295e-37 
        self.W = 1e7
        self.N_ph_k = (np.exp(((-1.60217656535e-22)/(self.temp * 1.38064852e-23))))**-1
        self.net_phonon_absorb = 0
    def initialize_system(self):
        for i in range(self.num_sites):
            self.site[i] = 1
        self.update_rate()
        for i in range(self.num_sites):
            if self.site[i] == 0:
                self.rate[i] = self.RATE_A2B
            else:
                self.rate[i] = self.RATE_B2A
        
            
    def compute_rate_sum(self):
        self.rate_sum = 0
        for i in range(self.num_sites):
            self.rate_sum += self.rate[i]
            
    def choose_site(self):
        self.rand_rate = random()
        while self.rand_rate == 0:
            self.rand_rate = random()
        x = 0
        rate_cum = self.rate[x]
        while abs(rate_cum) <= abs(self.rate_sum * self.rand_rate):
            rate_cum += (self.rate[int(x+1)])
            x += 1
        self.selected_site = int(x)
        
    def update_rate(self):
        zero_count = (self.site == 1).sum()
        k_count = (self.site == 0).sum()
        
        W_k_zero = self.W*(1 + self.N_ph_k)
        W_zero_k = self.W*(self.N_ph_k)
        
        self.RATE_B2A =  (-k_count * W_k_zero * (1 + zero_count)) + ((1 + k_count) * W_zero_k * zero_count) +(-self.gamma_k * k_count)
        self.RATE_A2B = -self.RATE_B2A + (-self.gamma_k * k_count)
        self.rate_listA2B.append(self.RATE_A2B)
        self.rate_listB2A.append(self.RATE_B2A)
        
    def execute_transition(self):
        if self.site[self.selected_site] == 0:
            self.site[self.selected_site] = 1
            self.net_phonon_absorb += -1
        elif self.site[self.selected_site] == 1:
            self.site[self.selected_site] = 0
            if (-self.gamma_k * ((self.site == 0).sum()))/self.RATE_B2A < random(): 
                self.net_phonon_absorb += 1
        self.update_rate()
        self.rate[self.selected_site] = (self.site[self.selected_site] == 1) and self.RATE_B2A or self.RATE_A2B
        
    def compute_time_increment(self):
        dt = 0.
        self.rand_time = random()
        dt = -np.log(self.rand_time) / abs(self.rate_sum)
        self.dt = dt
            
    def compute_pct_B(self):
        num_B = 0
        for i in range(self.num_sites):
            if self.site[i] == 1:
                num_B += 1
        self.pct_B = 100.0 * (num_B / self.num_sites)
        
    def evolve(self):
        i = 0
        for self.mcs in range(self.num_steps):
            self.compute_pct_B()
            self.pct_B_avg += self.pct_B
            self.B_list.append(self.pct_B)
            self.compute_rate_sum()
            self.choose_site()
            self.execute_transition()
            self.compute_time_increment()
            self.t += self.dt
            i += 1
            #print("This is Simulation Step ",i," of ", self.num_steps)
            
        self.pct_B_avg = self.pct_B_avg / self.num_steps
 
        
    def final_report(self):
        print("Simulation Summary:")
        print("Average Concentration of A States = ", 100.0 - self.pct_B_avg)
        print("Average Concentration of B States = ", self.pct_B_avg)
        print("Average Ratio of A to B States = ", (100.0 - self.pct_B_avg)/ self.pct_B_avg)
        print("Final Ratio of A to B Production Rates = ", self.RATE_B2A/self.RATE_A2B)
        print("# ---\n")
        print("Final Configuration \n")
        print(self.site)
        
        
        
def main():
    phonons = []
    for i in range(50):
        kmc = KMC(NUM_SITES, NUM_STEPS)
        kmc.initialize_system()
        kmc.evolve()
        #kmc.final_report()
        '''
        step = np.arange(kmc.mcs)
        step = (step / kmc.mcs) * kmc.t
        fig, ax = plt.subplots()
        #kmc.B_list = [100 - i for i in kmc.B_list]
        #ax.plot(step, kmc.rate_listA2B[0:kmc.mcs], label = r'A $\to$ B Rate')
        #ax.plot(step, kmc.rate_listB2A[0:kmc.mcs], label = r'B $\to$ A Rate')
        ax.plot(step, kmc.B_list[0:kmc.mcs], label = r'Ratio of $\frac{Lower State}{Higher State}$')
        ax.legend()
        ax.set_xlabel('Time (Seconds)')
        ax.set_ylabel(r'Percentage of Low Energy State ($\%$)')
        ax.set_title('Variable Rate Binary Kinetic Monte Carlo')
        plt.show()
        print("Total Phonons Absorbed, ", kmc.net_phonon_absorb)
        '''
        phonons.append(kmc.net_phonon_absorb)
        print("Iteration ",i," Completed")
    #fig, ax = plt.subplots()
    #ax.plot(np.arange(10), phonons, label = r'Phonons Absorbed Per Simulation')
    #ax.legend()
    #ax.set_xlabel('Simulation')
    #ax.set_ylabel('Number of Phonons Absorbed')
    #ax.set_title('Consistancy Testing')
    print(np.mean(phonons)," Phonons Absorbed on Average")
    #plt.show()
if __name__ == "__main__":
    main()
    
