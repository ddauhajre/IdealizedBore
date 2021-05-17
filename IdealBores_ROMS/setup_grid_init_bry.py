'''
Create ROMS grid, initial and boundary conditions
for idealized bore/SMC simulation

call the script like >> python setup_grid_init_bry.py params_file.py

where params_file.py is a separate script that initializes
parameters for grid, initial, and boundary conditions
'''
########################################
import os
import sys
from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib.pyplot as plt
import ROMS_depths as RD
from datetime import date
plt.ion()
######################################

#param_file = sys.argv[1]
#execfile('./'+param_file)
execfile('./params.py')



################
#Horizontal grid
################
x = np.arange(0,Lx,dx)
y = np.arange(0,Ly,dy)
nx = len(x)
ny = len(y)
pm = np.zeros([ny,nx])
pn = np.zeros([ny,nx])
angle = np.zeros([ny,nx]) #no rotation
mask_rho = np.zeros([ny,nx]) #all water
mask_rho[:,:] = 1.
pm[:,:] = 1./dx
pn[:,:] = 1./dy
x_rho, y_rho = np.meshgrid(x,y)
############################
#Create idealized bathmetry
#############################
h = np.zeros([ny,nx])
i_Ld = abs(x-Ld).argmin()
for j in range(ny):
   h[j,0:i_Ld] = H
   h[j,i_Ld::] = H - s_h*(x[i_Ld::]-Ld) 
h[h<hmin] = hmin
##############################################

#########################################
# Write and save variables to netcdf file
#########################################
grd_nc = netcdf(sim_name + '_grd.nc', 'w')
print ''
print 'Creating ROMS grid: ' + sim_name + '_grd.nc'
print ''

# SET GLOBAL ATTRIBUTES
grd_nc.title = 'ROMS grid produced by setup_grid_init_bry.py'
grd_nc.date = date.today().strftime("%B %d, %Y")

#SET DIMENSIONS
grd_nc.createDimension('xi_rho', nx)
grd_nc.createDimension('xi_u', nx-1)
grd_nc.createDimension('eta_rho', ny)
grd_nc.createDimension('eta_v', ny-1)

#CREATE AND WRITE VARIABLES

spherical = grd_nc.createVariable('spherical', 'c')
setattr(spherical, 'long_name', 'Grid type logical switch')
spherical[0]='F'



#xl_nc = grd_nc.createVariable('xl','f4', ('eta_rho', 'xi_rho'))
xl_nc = grd_nc.createVariable('xl','d')
setattr(xl_nc, 'long_name', 'domain length xi-direction')
setattr(xl_nc, 'units', 'meter')
xl_nc[0]=Lx


el_nc = grd_nc.createVariable('el','d')
setattr(el_nc, 'long_name', 'domain length in eta-direction')
setattr(el_nc, 'units', 'm')
el_nc[0]=Ly


xr_nc = grd_nc.createVariable('x_rho','f4', ('eta_rho', 'xi_rho'))
setattr(xr_nc, 'long_name', 'physical dimension in xi-direction')
setattr(xr_nc, 'units', 'm')
xr_nc[:,:] = x_rho


yr_nc = grd_nc.createVariable('y_rho','f4', ('eta_rho', 'xi_rho'))
setattr(yr_nc, 'long_name', 'physical dimension in eta-direction')
setattr(yr_nc, 'units', 'm')
yr_nc[:,:] = y_rho

h_nc = grd_nc.createVariable('h','f4', ('eta_rho', 'xi_rho'))
setattr(h_nc, 'long_name', 'bottom topography')
setattr(h_nc, 'units', 'm')
h_nc[:,:] = h

pm_nc = grd_nc.createVariable('pm','f4', ('eta_rho', 'xi_rho'))
setattr(pm_nc, 'long_name', 'curvilinear coordinate metric in XI-direction')
setattr(pm_nc, 'units', '1/m')
pm_nc[:,:] = pm

pn_nc = grd_nc.createVariable('pn','f4', ('eta_rho', 'xi_rho'))
setattr(pn_nc, 'long_name', 'curvilinear coordinate metric in ETA-direction')
setattr(pn_nc, 'units', '1/m')
pn_nc[:,:] = pn

f_nc = grd_nc.createVariable('f','f4', ('eta_rho', 'xi_rho'))
setattr(f_nc, 'long_name', 'Coriolis parameter at RHO-points')
setattr(f_nc, 'units', '1/s')
f_nc[:,:] = f0

angle_nc = grd_nc.createVariable('angle','f4', ('eta_rho', 'xi_rho'))
setattr(angle_nc, 'long_name', 'angle between east and XI-directions')
setattr(angle_nc, 'units', 'degrees')
angle_nc[:,:] = grid_ang

mask_nc = grd_nc.createVariable('mask_rho','f4', ('eta_rho', 'xi_rho'))
setattr(mask_nc, 'long_name', 'land/sea mask at rho-points')
setattr(mask_nc, 'units', 'land/water (0/1)')
mask_nc[:,:] = mask_rho

grd_nc.close()

print 'Saved grid file'


################################################
#Initial Condition
###############################################

print ' '
print '         Computing Initial Condition '

###################
#Make vertical grid
'''
Create approximate sigma level 
grid to mimic ROMS functionality
'''
####################
zeta = np.zeros([ny,nx])
Vtrans = 2
Vstret = 4
z_r = RD.set_depth(Vtrans,Vstret,theta_s,theta_b,hc,N,1,h,zeta)
z_w = RD.set_depth(Vtrans,Vstret,theta_s,theta_b,hc,N,5,h,zeta)
s,Cs_r = RD.stretching(Vstret,theta_s, theta_b, hc,N,0)
s,Cs_w = RD.stretching(Vstret,theta_s, theta_b, hc,N,1)
#Cs_r = Cs_r_temp[0:N]
####################################
#Idealized initial buoyancy profile
####################################
h_sbl = np.zeros([ny,nx])
for i in range(nx):
    h_sbl[:,i] = h0 + dh * np.exp( -((x[i] - x_fil) / L_fil)**2)


b = np.zeros([N,ny,nx])
for i in range(nx):
    for j in range(ny):
        a1 = Nb * (z_r[j,i,:] + H) 
        a2 = (0.5 * N0) * ( (1+ B) * z_r[j,i,:] - ( 1- B) *( h_sbl[j,i] + lambda_inv * np.log(np.cosh((1./lambda_inv) *(z_r[j,i,:] + h_sbl[j,i])))))
        a3 = (0.5 * N2) * ( (1+ B) * z_r[j,i,:] - ( 1- B) *( h2 + lambda_inv2 * np.log(np.cosh((1./lambda_inv2) *(z_r[j,i,:] + h2)))))
        b[:,j,i] = b0 + a1 + a2 - a3 


######################################
#Initialize temperature, not buoyancy
'''
T = b / (alpha * g)
where alpha = thermal expansion coefficient
'''
####################################
temp = b / (2e-4 * 9.81)
#####################
#Flow initilization
####################
u = np.zeros([N,ny,nx-1])
v = np.zeros([N,ny-1,nx])
ubar = np.zeros([ny,nx-1])
vbar = np.zeros([ny-1,nx])
zeta = np.zeros([ny,nx])


##########################
#Write Variables
#########################
init_nc = netcdf(sim_name + '_init.nc', 'w')
print ''
print 'Creating initial condition: ' + sim_name + '_init.nc'
print ''


# SET GLOBAL ATTRIBUTES
init_nc.title = 'ROMS initial conditon produced by setup_grid_init_bry.py'
init_nc.date = date.today().strftime("%B %d, %Y")
init_nc.VertCoordType='SM09'
init_nc.theta_s = theta_s
init_nc.theta_b = theta_b
init_nc.hc = hc
init_nc.Cs_r = Cs_r
init_nc.Cs_w = Cs_w

#SET DIMENSIONS
init_nc.createDimension('xi_rho', nx)
init_nc.createDimension('xi_u', nx-1)
init_nc.createDimension('eta_rho', ny)
init_nc.createDimension('eta_v', ny-1)
init_nc.createDimension('s_rho', N)
init_nc.createDimension('s_w', N+1)



#CREATE AND WRITE VARIABLES
otime_nc = init_nc.createVariable('ocean_time', 'd')
setattr(otime_nc, 'long_name', 'averaged time since initialization')
setattr(otime_nc, 'units', 'second')
otime_nc[0] = 0.

temp_nc = init_nc.createVariable('temp','f4', ('s_rho','eta_rho', 'xi_rho'))
setattr(temp_nc, 'long_name', 'potential temperature')
setattr(temp_nc, 'units', 'deg C')
temp_nc[:,:,:] = temp

salt_nc = init_nc.createVariable('salt','f4', ('s_rho','eta_rho', 'xi_rho'))
setattr(salt_nc, 'long_name', 'salinity')
setattr(salt_nc, 'units', 'psu')
salt_nc[:,:,:] = S0


u_nc = init_nc.createVariable('u','f4', ('s_rho','eta_rho', 'xi_u'))
setattr(u_nc, 'long_name', 'XI-velocity component')
setattr(u_nc, 'units', 'm/s')
u_nc[:,:,:] = u

ubar_nc = init_nc.createVariable('ubar','f4', ('eta_rho', 'xi_u'))
setattr(ubar_nc, 'long_name', 'barotropic XI-velocity component')
setattr(ubar_nc, 'units', 'm/s')
ubar_nc[:,:] = ubar

v_nc = init_nc.createVariable('v','f4', ('s_rho','eta_v', 'xi_rho'))
setattr(v_nc, 'long_name', 'ETA-velocity component')
setattr(v_nc, 'units', 'm/s')
v_nc[:,:,:] = v

vbar_nc = init_nc.createVariable('vbar','f4', ('eta_v', 'xi_rho'))
setattr(vbar_nc, 'long_name', 'barotropic ETA-velocity component')
setattr(vbar_nc, 'units', 'm/s')
vbar_nc[:,:] = vbar

zeta_nc = init_nc.createVariable('zeta','f4', ('eta_rho', 'xi_rho'))
setattr(zeta_nc, 'long_name', 'free-surface elevation')
setattr(zeta_nc, 'units', 'm')
zeta_nc[:,:] = zeta


init_nc.close()
###############################################

#################################
# SET UP BOUNDARY CONDITION 
'''
Internal tide wavemaker
that is function of
1st baroclinic mode of
western boundary stratification
'''
#################################
print ' '
print '         Computing Boundary Conditions '

import dynmodes as DM
tvec = np.arange(0,t_end, dt_bc)
nt = len(tvec)
u_west = np.zeros([nt,N,ny]) #inflow boundary condition of internal tide
ubar_west = np.zeros([nt,ny])
v_west = np.zeros([nt,N,ny-1])
vbar_west = np.zeros([nt,ny-1])
zeta_west = np.zeros([nt,ny])
temp_west = np.zeros([nt,N,ny])
salt_west = np.zeros([nt,N,ny])

dz_eig = 0.1 #vertical resolution used to solve TG eqn

#Only do modal decomposition at 1-y point
# as N^2=same everywhere
j = 0
Nsq = np.zeros([N+1])
Nsq[1:-1] = (b[1:,j,0] - b[:-1,j,0]) / (z_r[j,0,1:] - z_r[j,0,:-1])
Nsq[0] = Nsq[1]
Nsq[-1] = Nsq[-2]
z = z_w[j,0,:]
#Interpolate values to uniform grid for modal decomposition
z0 = z[0]
z1 = z[-1]
z_new = np.arange(z0,z1,dz_eig)
Nsq_in = np.interp(z_new, z, Nsq)
#Compute vertical modes of no-flow Taylor Goldstein system
wmodes, pmodes, rmodes, ce, dz, d2dz2 = DM.dynmodes(Nsq_in,z_new,1)

#Create U0(t) function to ramp up amplitude of tidal forcin
if U0_c_fixed:
   U0 = ce[0] * U0_c
#Create U0(t) function to ramp up amplitude of tidal forcing
U0_t = U0 * (tvec / t_full_U0)
U0_t[U0_t>U0] = U0


#pmodes_rho = 0.5 * (pmodes[:,1:] + pmodes[:,:-1])
pmodes_rho = np.interp(z_r[j,0,:],z_new,pmodes[0,:])
#Normalize pmodes_rho by its maximum
pmodes_rho = pmodes_rho / np.max(abs(pmodes_rho))
for j in range(ny):
    for k in range(N):
        u_west[:,k,j] = U0_t * np.sin(omega*tvec) * pmodes_rho[k]   

print '     U0 / c = ' + str(U0 / ce[0])

#Other fields match initial condition at offshore point
for t in range(nt):
    ubar_west[t,:] = ubar[:,0]
    v_west[t,:,:] = v[:,:,0]
    vbar_west[t,:] = vbar[:,0]
    zeta_west[t,:] = zeta[:,0]
    temp_west[t,:,:] = temp[:,:,0]
    salt_west[t,:,:] = S0 


##########################
#Write Variables
#########################
bry_nc = netcdf(sim_name + '_bry.nc', 'w')
print ''
print 'Creating boundary condition: ' + sim_name + '_bry.nc'
print ''


# SET GLOBAL ATTRIBUTES
bry_nc.title = 'ROMS boundary conditon produced by setup_grid_init_bry.py'
bry_nc.date = date.today().strftime("%B %d, %Y")
bry_nc.VertCoordType='SM09'
bry_nc.theta_s = theta_s
bry_nc.theta_b = theta_b
bry_nc.hc = hc
bry_nc.Cs_r = Cs_r
bry_nc.Cs_w = Cs_w

#SET DIMENSIONS
bry_nc.createDimension('xi_rho', nx)
bry_nc.createDimension('xi_u', nx-1)
bry_nc.createDimension('eta_rho', ny)
bry_nc.createDimension('eta_v', ny-1)
bry_nc.createDimension('s_rho', N)
bry_nc.createDimension('rec_time', None)


time_nc= bry_nc.createVariable('bry_time','d', ('rec_time' ))
setattr(time_nc, 'units', 'days')
time_nc[:] = tvec / 86400.

zeta_west_nc= bry_nc.createVariable('zeta_west','f4', ('rec_time','eta_rho' ), fill_value=1e33)
setattr(zeta_west_nc, 'units', 'deg C')
zeta_west_nc[:,:] = zeta_west

temp_west_nc= bry_nc.createVariable('temp_west','f4', ('rec_time', 's_rho','eta_rho' ), fill_value=1e33)
setattr(temp_west_nc, 'units', 'deg C')
temp_west_nc[:,:,:] = temp_west

salt_west_nc= bry_nc.createVariable('salt_west','f4', ('rec_time', 's_rho','eta_rho' ), fill_value=1e33)
setattr(salt_west_nc, 'units', 'psu')
salt_west_nc[:,:,:] = salt_west


u_west_nc= bry_nc.createVariable('u_west','f4', ('rec_time', 's_rho','eta_rho' ), fill_value=1e33)
setattr(u_west_nc, 'units', 'm/s')
u_west_nc[:,:,:] = u_west

ubar_west_nc= bry_nc.createVariable('ubar_west','f4', ('rec_time', 'eta_rho' ),fill_value=1e33)
setattr(ubar_west_nc, 'units', 'm/s')
ubar_west_nc[:,:] = ubar_west

v_west_nc= bry_nc.createVariable('v_west','f4', ('rec_time', 's_rho','eta_v' ),fill_value=1e33)
setattr(v_west_nc, 'units', 'm/s')
v_west_nc[:,:,:] = v_west

vbar_west_nc= bry_nc.createVariable('vbar_west','f4', ('rec_time', 'eta_v' ),fill_value=1e33)
setattr(vbar_west_nc, 'units', 'm/s')
vbar_west_nc[:,:] = vbar_west


bry_nc.close()

