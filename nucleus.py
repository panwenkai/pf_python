import numpy as np
def nuc(Nx, Ny, seed):
    phi = np.zeros((Nx, Ny))
    tempr = np.zeros((Nx, Ny))
    
    for i in range(Nx):
        for j in range(Ny):
            if ((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed):
                phi[i, j] = 1.0
                
    return phi, tempr

def nuc_arr(Nx, Ny, seed):
    phi = np.zeros((Nx+2, Ny+2))
    tempr = np.zeros((Nx+2, Ny+2))
    
    for i in range(Nx):
        for j in range(Ny):
            if ((i+1-Nx/2)*(i+1-Nx/2)+(j+1-Ny/2)*(j+1-Ny/2) < seed):
                phi[i+1, j+1] = 1.0
                
    return phi, tempr

def nuc_arr_eq(Nx, Ny, seed):
    phi = np.zeros((Nx+2, Ny+2))
    tempr = np.ones((Nx+2, Ny+2))
    
    for i in range(Nx):
        for j in range(Ny):
            if ((i+1-Nx/2) > (j+1-Ny/2)):
                phi[i+1, j+1] = 1.0
                
    return phi, tempr