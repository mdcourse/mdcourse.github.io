import random
import numpy as np

def initialization(n_atoms, T, dt, Lx, Ly, Lz):
    
    xlo, xhi = -Lx/2, Lx/2
    ylo, yhi = -Ly/2, Ly/2
    zlo, zhi = -Lz/2, Lz/2  

    # initialise positions and velocities
    x, y, z = np.zeros((3,n_atoms))
    xm, ym, zm = np.zeros((3,n_atoms))
    vx, vy, vz = np.zeros((3,n_atoms)) 

    for N in range(n_atoms):
        x[N] = random.random()*Lx+xlo
        y[N] = random.random()*Ly+ylo
        z[N] = random.random()*Lz+zlo
    
    vcom = np.zeros(3)
    ekin = 0
    for N in range(n_atoms):
        vx[N] = random.random()-0.5
        vy[N] = random.random()-0.5
        vz[N] = random.random()-0.5
        vcom += vx[N], vy[N], vz[N]
        ekin += vx[N]**2 + vy[N]**2 + vz[N]**2

    vcom /= n_atoms
    ekin /= n_atoms
    sfac = np.sqrt(3*T/ekin)

    for N in range(n_atoms):
        vx[N] = (vx[N]-vcom[0])*sfac
        vy[N] = (vy[N]-vcom[1])*sfac
        vz[N] = (vz[N]-vcom[2])*sfac
        xm[N] = x[N] - vx[N]*dt
        ym[N] = y[N] - vy[N]*dt
        zm[N] = z[N] - vz[N]*dt

    return x, y, z, vx, vy, vz, xm, ym, zm