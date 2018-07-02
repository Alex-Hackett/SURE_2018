# -*- coding: utf-8 -*-
"""
This is a translation of the example KMC C code found in 
doi:10.1016/j.cma.2008.03.010 into python 3.6
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

NUM_SITES = 24
NUM_STEPS = 64

RATE_A2B = 1.0
RATE_B2A = 5.0





class KMC:
    def __init__(self, sites_num, steps_num):
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
    def initialize_system(self):
        for i in range(self.num_sites):
            self.site[i] = 0
        for i in range(self.num_sites):
            self.prod[i] = 1
        for i in range(self.num_sites):
            self.rate[i] = RATE_A2B
            
    def compute_rate_sum(self):
        self.rate_sum = 0
        for i in range(self.num_sites):
            self.rate_sum += self.rate[i]
            
    def choose_site(self):
        self.rand_rate = random()
        x = 0
        rate_cum = self.rate[x]
        while rate_cum <= (self.rate_sum * self.rand_rate):
            rate_cum += self.rate[int(x+1)]
            x += 1
        self.selected_site = int(x)
        
    def execute_transition(self):
        site_copy = self.site[self.selected_site]
        self.site[self.selected_site] = self.prod[self.selected_site]
        self.prod[self.selected_site] = site_copy
        self.rate[self.selected_site] = (self.prod[self.selected_site] == 0) and RATE_B2A or RATE_A2B
        
    def compute_time_increment(self):
        dt = 0.
        self.rand_time = random()
        dt = -np.log(self.rand_time) / self.rate_sum
        self.dt = dt
            
    def compute_pct_B(self):
        num_B = 0
        for i in range(self.num_sites):
            if self.site[i] == 1:
                num_B += 1
        self.pct_B = 100.0 * (num_B / self.num_sites)
        
    def evolve(self):
        for self.mcs in range(self.num_steps):
            self.compute_pct_B()
            self.pct_B_avg += self.pct_B
            self.compute_rate_sum()
            self.choose_site()
            self.execute_transition()
            self.compute_time_increment()
            self.t += self.dt
            
        self.pct_B_avg = self.pct_B_avg / self.num_steps
        
    def final_report(self):
        print("Simulation Summary:")
        print("Average Concentration of A States = ", 100.0 - self.pct_B_avg)
        print("Average Concentration of B States = ", self.pct_B_avg)
        print("Average Ratio of A to B States = ", (100.0 - self.pct_B_avg)/ self.pct_B_avg)
        print("Ratio of A to B Production Rates = ", RATE_B2A/RATE_A2B)
        print("# ---\n")
        print("Final Configuration \n")
        print(self.site)
        
        
        
def main():
    kmc = KMC(NUM_SITES, NUM_STEPS)
    kmc.initialize_system()
    kmc.evolve()
    kmc.final_report()
    
if __name__ == "__main__":
    main()
    