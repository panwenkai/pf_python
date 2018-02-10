from nucleus import nuc_arr, nuc_arr_eq
import numpy as np
import matplotlib.pyplot as plt
import os

Nx = 300
Ny = 300;
NxNy = Nx*Ny

dx = 0.03
dy = 0.03
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
seed = 5.0

pix = 4.0*np.arctan(1.0)

top = 130
bottom = 170

# ---Initialize and introduce initial nuclei:
phi, tempr = nuc_arr_eq(Nx, Ny, seed)
lap_phi = np.zeros((Nx+2, Ny+2))
lap_tempr = np.zeros((Nx+2, Ny+2))
phidx = np.zeros((Nx+2, Ny+2))
phidy = np.zeros((Nx+2, Ny+2))
epsilon = np.zeros((Nx+2, Ny+2))
epsilon_deriv = np.zeros((Nx+2, Ny+2))

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
    
    
    phi[top,:] = phi[top+1,:]
    phi[bottom+1,:] = phi[bottom,:]
    phiold = np.array(phi)
    
    lap_phi[1:-1,1:-1] = (phi[0:-2,1:-1]+phi[2:,1:-1]+phi[1:-1,0:-2]+phi[1:-1,2:]-4.0*phi[1:-1,1:-1])/(dx*dy)
    lap_tempr[1:-1,1:-1] = (tempr[0:-2,1:-1]+tempr[2:,1:-1]+tempr[1:-1,0:-2]+tempr[1:-1,2:]-4.0*tempr[1:-1,1:-1])/(dx*dy)
    phidx[1:-1,1:-1] = (phi[2:,1:-1] - phi[0:-2,1:-1])/dx
    phidy[1:-1,1:-1] = (phi[1:-1,2:] - phi[1:-1,0:-2])/dy
    
    theta = np.arctan2(phidy, phidx)
    epsilon = epsilonb*(1.0+delta*np.cos(aniso*(theta-theta0)))
    epsilon_deriv = -epsilonb*aniso*delta*np.sin(aniso*(theta-theta0))
    
    
    term1 = (epsilon[1:-1,2:]*epsilon_deriv[1:-1,2:]*phidx[1:-1,2:] - epsilon[1:-1,0:-2]*epsilon_deriv[1:-1,0:-2]*phidx[1:-1,0:-2])/dy
    term2 = -(epsilon[2:,1:-1]*epsilon_deriv[2:,1:-1]*phidy[2:,1:-1] - epsilon[0:-2,1:-1]*epsilon_deriv[0:-2,1:-1]*phidy[0:-2,1:-1])/dx        
    m = alpha/np.pi * np.arctan(gamma*(teq-tempr[1:-1,1:-1]))
    phi[1:-1,1:-1] = phi[1:-1,1:-1] +(dtime/tau)*(term1+term2+epsilon[1:-1,1:-1]**2*lap_phi[1:-1,1:-1] + phiold[1:-1,1:-1]*(1.0-phiold[1:-1,1:-1])*(phiold[1:-1,1:-1] -0.5 + m))
    phi[top,:] = phi[top+1,:]
    phi[bottom+1,:] = phi[bottom,:]
    phi[:top,:] = phiold[:top,:]
    phi[bottom+2:,:] = phiold[bottom+2:,:]
    tempr[1:-1,1:-1] = tempr[1:-1,1:-1] +dtime*lap_tempr[1:-1,1:-1] + kappa*(phi[1:-1,1:-1]-phiold[1:-1,1:-1])
    
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