#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
from random import random

Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.001

dx = Lx/Nx

N = 20

print Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
tab = np.exp(-SpaceX**2)
tab /= np.sqrt(np.linalg.norm(tab))
tab += 1

tab = tab.astype(np.complex)
tab2 = tab.copy()

Kx = np.arange(Nx) - np.append(np.zeros(Nx/2),np.ones(Nx/2)*Nx)
Kx *= 2*np.pi/Lx

K2 = np.exp(-1j*Kx**2*dt)

Tab = []
Tab2 = []

E = np.ones(Nx)

for i in range(6):
    tab = ifft(K2*fft(tab))
    r = np.abs(tab2)
    #source of the problem
    E = tab2/np.abs(tab2)
    de = np.angle(E/np.roll(E, 1))/dx
    d2e = np.angle(E/np.roll(E, 1)**2*np.roll(E, 2))/dx/dx

    dr = (r - np.roll(r, 1))/dx
    d2r = (np.roll(r, -1) - 2*r + np.roll(r, 1))/dx/dx
    dPsi = 1j*d2r * E - tab2*d2e - 2*dr*de * E - 1j*tab2*de**2
    dPsi *= dt
    tab2 += dPsi
    Tab.append(np.abs(tab))
    Tab2.append(tab2)
    print i, np.linalg.norm(tab - tab2)


plt.plot(fft(de))
plt.plot(fft(d2e))
plt.show()

exit(0)

Tab = np.array(Tab)
plt.imshow(Tab.T, cmap='gnuplot')
plt.show()