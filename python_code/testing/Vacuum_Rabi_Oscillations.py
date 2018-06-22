#!/usr/bin/python
'''
quick implementation of the J-C model of vacuum oscillations based
on the code provided in the QuTiP documentation
'''

#Imports
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import astropy as ap
#Defining parameters
#hbar = 1 units used

#Simulation time (picoseconds)
time = 100
#Simulation steps
inter = 100000

#Cavity Mode Freq
wc = 1.0 * 2 * np.pi
#Atom Level Freq
wa = 1.0 * 2 * np.pi
#Coupling constant
g = 0.05 * 2 * np.pi
#Cavity dissipation rate
kappa = 0.0#0.005
#atomic dissipation rate
gamma = 0.0#0.05
#Number of cavity mode Fock States
N = 15
#Temperature, in freq units
n_th_a = 0.0
#Use the rotating wave approximation
use_rwa = True

tlist = np.linspace(0,time,inter)

#Generate the initial excited atomic Fock state

psi0 = tensor(basis(N,0), basis(2,1))

#Generate Operators

#Photon field annihilation operator
a = tensor(destroy(N), qeye(2))
#matter field annihilation operator
sm = tensor(qeye(N), destroy(2))

#Setup the Hamiltonian
if use_rwa:
    H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
else:
    H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag())


#Operator set required to describe dissipations

#create a list to hold the collapse operators
c_op_list = []
#And describe the dissipation rate of the states
rate = kappa * (1 + n_th_a)
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * a)

rate = kappa * n_th_a
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * a_dag())

rate = gamma
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * sm)


'''
Time evolve the system using the Lindblad equation. Get the expectation values of the photonic and matter field number operators, this can just be passed as an argument to the built in equation solver.
'''

numb_op = [a.dag() * a, sm.dag() * sm]
output = mesolve(H, psi0, tlist, c_op_list, numb_op)

'''
Visualization
Plot the  excitation probabilities of the cavity and the atomic system as a function of time
'''
'''
fig, ax = plt.subplots()
ax.plot(tlist, output.expect[0], label = 'Cavity')
ax.plot(tlist, output.expect[1], label = 'Atom Excited State')
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Occupation Probability')
ax.set_title('Vacuum Rabi Oscillations')
plt.show()
'''

'''
Test of detuning the cavity mode coupling
'''










