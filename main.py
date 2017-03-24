#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.0001

dx = Lx/Nx

N = 500

print Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
tab = np.exp(-SpaceX**2)
tab /= np.sqrt(np.linalg.norm(tab))
tab += 1

tab = tab.astype(np.complex)
tab2 = tab.copy()

Kx = np.arange(Nx) - np.append(np.zeros(Nx/2), np.ones(Nx/2)*Nx)
Kx *= 2*np.pi/Lx

K2 = np.exp(-1j*Kx**2*dt)

Tab = []
Tab2 = []

E = np.ones(Nx)

for i in range(N):
    tab = ifft(K2*fft(tab))
    E = tab/np.abs(tab)
    dTheta = np.imag(np.roll(E, 1)/E - 1)/dx
    d2Theta = (dTheta - np.roll(dTheta, 1))/dx

    r = np.abs(tab2)
    dr = (r - np.roll(r, 1))/dx
    d2r = (dr - np.roll(dr, 1))/dx

    # tab2 = np.exp(-1j*dt*((dTheta)**2 - 1j*d2Theta))*tab2
    #i don't know why but without 1j thing it's more accurate
    # dPsi = -(d2r + 2j*dTheta*dr)*E
    #again, when I remove j it coverges better
    tab2 = np.exp(-1j*dt*((dTheta)**2 - d2Theta))*tab2
    dPsi = d2r*E*1j*dt
    tab2 += dPsi



    Tab.append(np.abs(tab))
    Tab2.append(tab2)
    print i, np.linalg.norm(tab - tab2)


exit(0)

Tab = np.array(Tab)
plt.imshow(Tab.T, cmap='gnuplot')
plt.show()
