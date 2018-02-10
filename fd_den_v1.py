from nucleus import nuc
import numpy as np
import matplotlib.pyplot as plt

Nx = 300
Ny = 300;
NxNy = Nx*Ny

dx = 0.03
dy = 0.03
nstep = 4000
nprint= 50
dtime = 1.0e-4

# --- Material specific parameters:

tau = 0.0003
epsilonb = 0.01
mu = 1.0
kappa = 1.8
delta = 0.02
aniso = 4.0
alpha = 0.9
gamma = 10.0
teq = 1.0
theta0 = 0.2
seed = 5.0

pix = 4.0*np.arctan(1.0)

# ---Initialize and introduce initial nuclei:
phi, tempr = nuc(Nx, Ny, seed)
lap_phi = np.zeros((Nx, Ny))
lap_tempr = np.zeros((Nx, Ny))
phidx = np.zeros((Nx, Ny))
phidy = np.zeros((Nx, Ny))
epsilon = np.zeros((Nx, Ny))
epsilon_deriv = np.zeros((Nx, Ny))

phi_path = "D://phi//"
tempr_path = "D://tempr//"

for istep in range(nstep):
    print istep
    filename = str(istep)
    filename += ".jpg"
    
    for i in range(Nx):
        for j in range(Ny):
            
            jp = j + 1
            jm = j - 1
            
            ip = i + 1
            im = i - 1
            
            if (im == -1):
                im = Nx - 1
                
            if (ip == Nx):
                ip = 0
                
            if (jm == -1):
                jm = Ny - 1
               
            if (jp == Nx):
                jp = 0
                
            hne = phi[ip,j]
            hnw = phi[im,j]
            hns = phi[i,jm]
            hnn = phi[i,jp]
            hnc = phi[i,j]

            lap_phi[i,j] = (hnw + hne + hns + hnn -4.0*hnc)/(dx*dy)

            hne = tempr[ip,j]
            hnw = tempr[im,j]
            hns = tempr[i,jm]
            hnn = tempr[i,jp]
            hnc = tempr[i,j]

            lap_tempr[i,j] = (hnw + hne + hns + hnn -4.0*hnc)/(dx*dy)
            
            phidx[i,j] = (phi[ip,j] - phi[im,j])/dx
            phidy[i,j] = (phi[i,jp] - phi[i,jm])/dy
            
            theta = np.arctan2(phidy[i,j], phidx[i,j])
            
            epsilon[i,j] = epsilonb*(1.0+delta*np.cos(aniso*(theta-theta0)))
            epsilon_deriv[i,j] = -epsilonb*aniso*delta*np.sin(aniso*(theta-theta0))
            
            
    for i in range(Nx):
        for j in range(Ny):
            
            jp = j + 1
            jm = j - 1
            
            ip = i + 1
            im = i - 1
            
            if (im == -1):
                im = Nx - 1
                
            if (ip == Nx):
                ip = 0
                
            if (jm == -1):
                jm = Ny - 1
               
            if (jp == Nx):
                jp = 0
                
            phiold = phi[i,j]
            
            term1 = (epsilon[i,jp]*epsilon_deriv[i,jp]*phidx[i,jp] - epsilon[i,jm]*epsilon_deriv[i,jm]*phidx[i,jm])/dy
            term2 = -(epsilon[ip,j]*epsilon_deriv[ip,j]*phidy[ip,j] - epsilon[im,j]*epsilon_deriv[im,j]*phidy[im,j])/dx;
            
            m = alpha/np.pi * np.arctan(gamma*(teq-tempr[i,j]))
            
            phi[i,j] = phi[i,j] +(dtime/tau)*(term1+term2+epsilon[i,j]**2*lap_phi[i,j] + phiold*(1.0-phiold)*(phiold -0.5 + m))
            
            tempr[i,j] = tempr[i,j] +dtime*lap_tempr[i,j] + kappa*(phi[i,j]-phiold)
    fig = plt.figure(frameon=False)
    plt.imshow(phi)
    plt.colorbar()
    fig.savefig(phi_path + filename, dpi=500)
    fig.clf()
    plt.close()
    
    fig = plt.figure(frameon=False)
    plt.imshow(tempr)
    plt.colorbar()    
    fig.savefig(tempr_path + filename, dpi=500)
    fig.clf()
    plt.close()