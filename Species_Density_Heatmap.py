"""
We want to use differential equations to model
the number density of each species (e.g. )

https://arxiv.org/pdf/2609.01366
2.3.3 Secondary heating: dd example

"""

# importing packages
import numpy as np
import matplotlib.pyplot as plt
import plasmapy
import csv

# intializing variables

total_t = 300.0  # total time in seconds
dt = 0.01  # timestep in seconds


b_t = 0.5  # branching ratio for tritons
b_3 = 0.5  # branching ratio for 3He nuclei



# integrated cross-section
# x = (energy/energy threshold) and the values for a0, a1, a2, a3, a4 are the fits
# energy threshold and fits can be found in el-D2-X1Sg/excitation/fits/

sigma = ((x-1)/x * (a0^2/x + a1/x^2 + a2/x^3 + a3/x^4 + a4/x^5))

# D-D Reactivity Fit Parameters

sigmav_dd
sigmav_dt
sigmav_d3

R_dd = 0.5 * ((n_d)^2) * (sigmav_dd)

S_t = (b_t) * (R_dd)
S_t = (b_3) * (R_dd)