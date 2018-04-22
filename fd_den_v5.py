from nucleus import nuc_arr, nuc_arr_eq
import numpy as np
import matplotlib.pyplot as plt
import os

Nx = 300
Ny = 300;
NxNy = Nx*Ny

# Here we add more physical parameters
sigma = 0.4 #J/m^2
T_m = 1687 #K
mu = 0.07 #m/(K s)
L = 4.1474e-9 # J/m^3
delta = 1e-9 #m

dx = 1e-9
dy = 1e-9
nstep = 100000
nprint= 50
dtime = 1.0e-4

# --- Material specific parameters:

tau = 0.00015
epsilonb = 0.01
mu = 1.0
kappa = 1.8
delta = 0.02
aniso = 0.0
alpha = 0.45
gamma = 10.0
teq = 1.0
theta0 = 0.2
top = 130
bottom = 170
seed = 5.0
# ---Initialize and introduce initial nuclei:
phi, tempr = nuc_arr_eq(Nx, Ny, seed)
lap_phi = np.zeros((Nx, Ny))
lap_tempr = np.zeros((Nx, Ny))
phidx = np.zeros((Nx+2, Ny+2))
phidy = np.zeros((Nx+2, Ny+2))
epsilon = np.zeros((Nx, Ny))
epsilon_deriv = np.zeros((Nx, Ny))

# when phi is 0, we are in liquid
# when phi is 1, we are in solid
heat_capacity_0 = np.ones((Nx+2,Ny+2)) * 2.2e6 # J/ (m^3 K)
heat_capacity_1 = np.ones((Nx+2,Ny+2)) * 2.4e6 # J/ (m^3 K)

heat_conductivity_0 = 580 # W/ (m K)
heat_conductivity_1 = 200 # W/ (m K)

phi_path = "D://phi_6//"
tempr_path = "D://tempr_6//"

if not os.path.exists(phi_path):
    os.makedirs(phi_path)
if not os.path.exists(tempr_path):
    os.makedirs(tempr_path)

for istep in range(nstep):
    
    filename = str(istep)
    filename += ".jpg"
    
    phi[0,:] = phi[-2,:]
    phi[-1,:] = phi[1,:]
    phi[:,0] = phi[:,-2]
    phi[:,-1] = phi[:,1]
    
    tempr[0,:] = tempr[-2,:]
    tempr[-1,:] = tempr[1,:]
    tempr[:,0] = tempr[:,-2]
    tempr[:,-1] = tempr[:,1]
    phi[:top,:] = np.zeros((top,Ny+2))
    phi[bottom+2:,:] = np.zeros((Nx-bottom,Ny+2))
    
    heat_conductivity = heat_conductivity_0 + phi*heat_conductivity_1
    tempr_y_gredient = tempr[:-1,:] - tempr[1:,:]
    heatflow_y = 2*(heat_conductivity[:-1,:]*heat_conductivity[1:,:])*tempr_y_gredient/(heat_conductivity[:-1,:]+heat_conductivity[1:,:])
    tempr_x_gredient = tempr[:,:-1] - tempr[:,1:]
    
    phi[top,:] = phi[top+1,:]
    phi[bottom+1,:] = phi[bottom,:]
    phiold = np.array(phi)
    
    lap_phi = (phi[0:-2,1:-1]+phi[2:,1:-1]+phi[1:-1,0:-2]+phi[1:-1,2:]-4.0*phi[1:-1,1:-1])/(dx*dy)
    lap_tempr = (tempr[0:-2,1:-1]+tempr[2:,1:-1]+tempr[1:-1,0:-2]+tempr[1:-1,2:]-4.0*tempr[1:-1,1:-1])/(dx*dy)
    term1 = lap_phi - phi[1:-1,1:-1]*(1-phi[1:-1,1:-1])*(1-2*phi[1:-1,1:-1])/(delta**2)
    term1 *= sigma * T_m * mu /L
    term2 = (T_m - tempr[1:-1,1:-1])*phi[1:-1,1:-1]*(1-phi[1:-1,1:-1])
    term2 *= (-1*mu/delta)
    phi[1:-1,1:-1] = phi[1:-1,1:-1] + dtime*(term1+term2)
    phi[top,:] = phi[top+1,:]
    phi[bottom+1,:] = phi[bottom,:]
    phi[:top,:] = phiold[:top,:]
    phi[bottom+2:,:] = phiold[bottom+2:,:]
    tempr[1:-1,1:-1] = tempr[1:-1,1:-1] +dtime*lap_tempr + kappa*(phi[1:-1,1:-1]-phiold[1:-1,1:-1])
    
    if istep % 100 != 0:
        continue
    print istep
    fig = plt.figure(frameon=False)
    plt.imshow(phi[1:-1,1:-1])
    plt.colorbar()
    fig.savefig(phi_path + filename, dpi=500)
    fig.clf()
    plt.close()
    
    fig = plt.figure(frameon=False)
    plt.imshow(tempr[1:-1,1:-1])
    plt.colorbar()    
    fig.savefig(tempr_path + filename, dpi=500)
    fig.clf()
    plt.close()