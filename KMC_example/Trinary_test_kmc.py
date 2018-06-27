# -*- coding: utf-8 -*-
"""
This is an attempt to move from a binary state kinetic monte carlo code to an
arbitrary, Potts model style code, recalling that the quantum well confinment
in the actual physical system will result in the quantization of the possible
polariton wavevectors, thus this is an acceptable physical model
let's try three different states
Let's make all the rates just 1 for now
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

NUM_SITES = 1000
NUM_STEPS = 100000






class KMC:
    def __init__(self, sites_num, steps_num):
        self.num_sites = sites_num
        self.num_steps = steps_num
        self.site = np.zeros(sites_num, dtype = int)
        self.mcs = 0
        self.rate = np.zeros((self.num_sites,3), dtype = float)
        self.rate_sum = 0.
        self.t = 0.
        self.pct_B_avg = 0.
        self.pct_C_avg = 0.
        self.pct_A_avg = 0.
        self.rand_rate = 0.
        self.rand_time = 0.
        self.selected_site = 0
        self.pct_B = 0.
        self.pct_C = 0.
        self.dt = 0.
        self.B_list = []
        self.C_list = []
        self.A_list = []
        
    def initialize_system(self):
        '''
        num_sites by 3 matrix, first row is the A<-->B Rate, second row
        is the A<-->C rate, 3rd is B<-->C rate
        '''
        for i in range(self.num_sites):
            self.site[i] = 0
        for i in range(self.num_sites):
            for j in range(3):
                self.rate[i, j] = 1.0
            
        
            
    def compute_rate_sum(self):
        self.rate_sum = 0
        for i in range(self.num_sites):
            for j in range(3):
                self.rate_sum += self.rate[i,j]
            
    def choose_event(self):
        self.rand_rate = random()
        x = 0
        x_1 = 0
        x_2 = 0
        rate_cum = self.rate[x,x]
        for i in range(self.num_sites):
            for j in range(3):
                if rate_cum <= (self.rate_sum * self.rand_rate):
                    rate_cum += self.rate[i,j]
                    x_1 = i
                    x_2 = j
                else:
                    break
        self.selected_site = [x_1,x_2]
        
        
    def update_rate(self):
        self.RATE_B2A = 10*(np.sqrt(((self.site == 1).sum())) + 0.005)
        self.RATE_A2B = 10*(np.sqrt(((self.site == 0).sum())) + 0.005)
        self.rate_listA2B.append(self.RATE_A2B)
        self.rate_listB2A.append(self.RATE_B2A)
        
    def execute_transition(self):
        if self.selected_site[1] == 0:
            if self.site[self.selected_site[0]] == 0:
                self.site[self.selected_site[0]] = 1
            else:
                self.site[self.selected_site[0]] = 0
                
        elif self.selected_site[1] == 1:
            if self.site[self.selected_site[0]] == 0:
                self.site[self.selected_site[0]] = 2
            else:
                self.site[self.selected_site[0]] = 0
            
        elif self.selected_site[1] == 2:
            if self.site[self.selected_site[0]] == 1:
                self.site[self.selected_site[0]] = 2
            else:
                self.site[self.selected_site[0]] = 1
            
       
        
    def compute_time_increment(self):
        dt = 0.
        self.rand_time = random()
        dt = -np.log(self.rand_time) / self.rate_sum
        self.dt = dt
            
    def compute_pct(self):
        num_B = 0
        for i in range(self.num_sites):
            if self.site[i] == 1:
                num_B += 1
        self.pct_B = 100.0 * (num_B / self.num_sites)
        self.B_list.append(self.pct_B)
        num_C = 0
        for i in range(self.num_sites):
            if self.site[i] == 2:
                num_C += 1
        self.pct_C = 100.0 * (num_C / self.num_sites)
        self.C_list.append(self.pct_C)
        self.A_list.append(100 - self.pct_C - self.pct_B)
        
    def evolve(self):
        i = 0
        for self.mcs in range(self.num_steps):
            self.compute_pct()
            self.compute_rate_sum()
            self.choose_event()
            self.execute_transition()
            self.compute_time_increment()
            self.t += self.dt
            print("Simulation Time = ",self.t," Seconds")
        self.pct_B_avg = self.pct_B / self.num_steps
        self.pct_C_avg = self.pct_C / self.num_steps

 
        
    def final_report(self):
        print("Simulation Summary:")
        print("Average Concentration of A States = ", 100.0 - self.pct_B_avg)
        print("Average Concentration of B States = ", self.pct_B_avg)
        print("Average Ratio of A to B States = ", (100.0 - self.pct_B_avg)/ self.pct_B_avg)
        print("Final Ratio of A to B Production Rates = ", self.RATE_B2A/self.RATE_A2B)
        print("# ---\n")
        print("Final Configuration \n")
        print(self.site)
        
        
        

kmc = KMC(NUM_SITES, NUM_STEPS)
kmc.initialize_system()
kmc.evolve()
#kmc.final_report()

step = np.arange(kmc.mcs)
step = (step / kmc.mcs) * kmc.t
fig, ax = plt.subplots()
ax.plot(step, kmc.B_list[0:kmc.mcs], label = r'Percentage of B Sites')
ax.plot(step, kmc.C_list[0:kmc.mcs], label = r'Percentage of C Sites')
ax.plot(step, kmc.A_list[0:kmc.mcs], label = r'Percentage of A Sites')
ax.legend()
ax.set_xlabel('Time (Seconds)')
ax.set_ylabel(r'Percentage of Sites ($\%$)')
ax.set_title('Trinary Kinetic Monte Carlo')
plt.show()


    
