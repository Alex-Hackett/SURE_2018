# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:12:49 2018

@author: Alex
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

import polariton_scatter_pack as psp

hbar = 6.582119514e-16


