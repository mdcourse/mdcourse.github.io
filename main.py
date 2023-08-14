import numpy as np
import sys

n_atoms = 10
T = 1.0
dt = 0.005
box = np.array([20, 20, 20])
LJ_cutoff = 3
delta_dump = 10
Lx, Ly, Lz = 10, 10, 10

sys.path.append("functions")
from initialization import initialization

initialization(n_atoms, T, dt, Lx, Ly, Lz)