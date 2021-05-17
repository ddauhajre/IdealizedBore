'''
Declare all parameter names
and setup netcdf files for each set
'''
import numpy as np

base_name = 'VariableStrat' 
nsims =10
b0_l = [2.5e-2 for r in range(nsims)]
B_l = [0.025 for r in range(nsims)]
lambda_inv_l = [5. for r in range(nsims)]

Nb_s = 2e-4
Nb_w = 1e-6
N0_s =0
N0_w = 2e-4
N2_s =0
N2_w = 1.6e-4
ddb = (Nb_s - Nb_w) / (nsims)
dd0 = (N0_w - N0_s) / (nsims)
dd2 = (N2_w - N2_s) / (nsims)
Nb_l = np.arange(Nb_w,Nb_s+ddb,ddb)
N0_l = np.arange(N0_s,N0_w+dd0,dd0)[::-1]
N2_l = np.arange(N2_s,N2_w+dd2,dd2)[::-1]


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












