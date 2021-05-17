'''
Parameters for grid, initial and boundary conditions
'''
#sim_name = 'PycTest'
###########################
#Parameters
##########################
#Grid Paramters
dx = 35 #meters
dy = 35 #meters
Lx = 21140
Ly = 2310 
N = 150
theta_s = 6.
theta_b = 6.
hc = 60.

#Bathymetry Paramters
H = 170 #total tdepth
s_h = 8e-3 #bathymetric slope
Ld = 150 #length of flat shelf
hmin = 2. #minimum depth

#Other grid parameters
f0 = 8.3e-5 #coriolis parameter
#f0 = 0.
grid_ang = 0. #grid rotation (default=0-->xi=east)


x_fil = 10e3 #even if dh=0, need these to be >0 
L_fil = 1e3 #even if dh=0, need these >0
h0 = 50 
dh = 0 
S0 = 32 #constant salinity
########################
#FLOW INITILIZATION
'''
rest ==>, u,v,w,zeta =0 at t=0
'''
#########################
flow_init = 'rest'
#####################################
# Boundary Condition Parameters
#####################################
#Time of run
t_end = 10*24*3600. #seconds
dt_bc = 10*60 #time-step for boundary conditions
t_full_U0 = 24.*2 * 3600 #(seconds) time for full U0 forcing
omega = 1.41e-4 #forcing frequency

#Options to fix wave amplitude (U0 / c) based on modal speed (c)
U0_c_fixed = False
if U0_c_fixed:
   U0_c = 0.1
else:
   U0 = 0.08 #amplitude of internal tidal forcing

