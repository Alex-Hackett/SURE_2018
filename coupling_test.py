#!/usr/bin/python
'''
quick implementation of the J-C model of vacuum oscillations based
on the code provided in the QuTiP documentation
'''

#Imports
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from qutip import *

#Defining parameters
#hbar = 1 units used

def cavitySim(time, kappa, gamma, wc, wa):
    '''
    Simulations the time evolution of an atomic Fock state coupled to a 15
    photonic cavity modes, for a time "time" picoseconds, in 1/1000 ps 
    intervals. Kappa and Gamma are the photonic and atomic dissipation factors
    respectively and wc and wa are the photonic and atomic frequencies in 
    units of 2pi
    '''
    

    #Simulation time (picoseconds)
    #Simulation steps
    inter = time*1000
    
    #Cavity Mode Freq
    wc = wc * 2 * np.pi
    #Atom Level Freq
    wa = wa * 2 * np.pi
    #Coupling constant
    g = 0.05 * 2 * np.pi
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
    
    excited_cavity_occupation = output.expect[0]
    excited_matter_occupation = output.expect[1]
    return excited_cavity_occupation, excited_matter_occupation



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


detuning = np.linspace(0, 2, 100)
cavity_list = []
matter_list = []
for i in detuning:
    cavity, matter = cavitySim(2, 0.0, 0.0, 1.0, 1.0 + i)
    cavity = max(cavity)
    matter = max(matter)
    cavity_list.append(cavity)
    matter_list.append(matter)
    print('System Evolved for Detuning of ',i)
    
fig, ax = plt.subplots()
ax.plot(detuning, cavity_list, label = 'Cavity')
ax.plot(detuning, matter_list, label = 'Matter')
ax.legend()
ax.set_xlabel(r'$\Delta$, Detuning')
ax.set_ylabel('Maximum Excited State Occupation Probability')
ax.set_title('Effect of Detuning on Vacuum Rabi Oscillations')
plt.show()
    




