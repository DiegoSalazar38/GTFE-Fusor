# The plasma we are working with is a good approximation of a perfect gas (low density, low pressure high temperature)
# Goal is to plot the velocity (or kinetic energy) distribution function for the ions in the plasma
# We want to compare the Maxwell-Boltzmann distribution to our result

#importing our packages
import plasmapy
import numpy as np
import matplotlib.pyplot as plt

#Matplotlib syntax; creates our plots
fig, ax = plt.subplots()

num_points = 10**(3)
KE_MAX = 50 #arbitrary max for our plot, in eV
N_MAX = 1 #

#probability density defined by 
#p_density = np.arange(0,P_MAX)
k_energy = np.arange(0,KE_MAX, KE_MAX/num_points)

N = 10**(9) #number of particles
k = 8.617 * 10**(-5) #in eV/Kelvin
m = 1875.612 * 10**(6) #in eV/c^2
e = 1.602 * 10**(-19) #in Coulombs

v = (2*k_energy/m)**(1/2) #in m/s
P = 0.2666 * 1/e #pressure in eV/m^3
V = np.pi*0.25**2 #in m^3
T = P*V/N*k # Temperature using perfect gas approximation, in Kelvin

#From Maxwell-Boltzmann distribution equation
#N_density = 4*np.pi*N*(m/(2*np.pi*k*T))**(3/2)*v**(2)*np.exp(-m*v**2/(2*k*T))
N_density = 2*(k_energy/np.pi)**(0.5)*(k*T)**(-3/2)*np.exp(-k_energy/(k*T))
normalized_N_density = N_density / N #makes the integral of the distribution function over the entire domain equal to 1
#We're plotting probability density vs kinetic energy of our ions
ax.plot(k_energy, N_density)
plt.xlabel("Kinetic Energy (KE, eV)")
plt.ylabel("Fraction of Total Particles (N, eV^-1)")
plt.show()


