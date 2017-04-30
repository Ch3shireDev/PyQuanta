#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
from numpy import conj
from numpy.linalg import norm
import matplotlib.pyplot as plt

N = 500
Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.000001

dx = Lx/Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
tab = np.exp(-SpaceX**2 + 0.1j*SpaceX)
tab /= np.sqrt(np.linalg.norm(tab))

tab = tab.astype(np.complex)
tab2 = tab.copy()
tab3 = tab.copy()

Kx = 2j*np.pi/Lx*(np.arange(Nx) - np.append(np.zeros(Nx/2), np.ones(Nx/2)*Nx))

K2 = np.exp(-1j*Kx**2*dt)


for i in range(N):
    tab = ifft(K2*fft(tab))
    tab2 += -1j*dt*ifft(Kx**2*fft(tab2))
    print i, norm(tab-tab2)