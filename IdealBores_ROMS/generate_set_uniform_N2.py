'''
Declare all parameter names
and setup netcdf files for each set
'''
import numpy as np

base_name = 'UniformStrat' 
nsims =15
b0_l = [2.5e-2 for r in range(nsims)]
B_l = [0.025 for r in range(nsims)]
lambda_inv_l = [5. for r in range(nsims)]

Nb_l = [9.7e-7, 1.9e-6, 8.2e-6, 1.5e-5, 1.9e-5, 3.3e-5, 4e-5, 7.9e-5, 1.2e-4
        , 1.6e-4, 2e-4, 2.3e-4, 2.7e-4, 3.1e-4, 3.5e-4] 
N0_l = [0 for n in range(nsims)] 
N2_l = [0 for n in range(nsims)] 

h2_l = [70 for r in range(nsims)]
lambda_inv2_l = [5. for r in range(nsims)]


for n in range(nsims):
    sim_name = base_name + str(n+1)
    b0 = b0_l[n]
    B = B_l[n]
    lambda_inv = lambda_inv_l[n]
    Nb = Nb_l[n]
    N0 = N0_l[n]
    N2 = N2_l[n]
    h2 = h2_l[n]
    lambda_inv2 = lambda_inv2_l[n]
    execfile('./setup_grid_init_bry.py')












