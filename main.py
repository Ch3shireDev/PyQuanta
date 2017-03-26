#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.001

dx = Lx/Nx


print Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
tab = np.exp(-SpaceX**2)
tab /= np.sqrt(np.linalg.norm(tab))
tab += 0.9

tab = tab.astype(np.complex)
tab2 = tab.copy()

Kx = np.arange(Nx) - np.append(np.zeros(Nx/2), np.ones(Nx/2)*Nx)
Kx *= 2*np.pi/Lx

K2 = np.exp(1j*Kx**2*dt)

Tab = []
Tab2 = []




N = 200
E = np.ones(Nx).astype(np.complex)

r = np.abs(tab2)

for i in range(N):
    tab = ifft(K2*fft(tab))

    # E = tab/np.abs(tab)
    dTheta = np.real(1j*(1-np.roll(E, 1)/E))/dx
    d2Theta = (dTheta - np.roll(dTheta, 1))/dx

    dr = (r - np.roll(r, 1))/dx
    d2r = (dr - np.roll(dr, 1))/dx

    #i dPsi = (dTheta**2-id2Theta)*Psi - (d2r + 2i dTheta*dR)*E

    r = np.exp(-(2*dTheta*dr/r + d2Theta)*dt)*r
    # E = (r+1j*(dTheta**2*r)*dt)/np.sqrt(r**2 + dTheta**4*r**2)*E
    E = tab/np.abs(tab)

    tab2 = r*E

    Tab.append(np.angle(E))

    print i, np.linalg.norm(tab - tab2)

exit(0)
Tab = np.array(Tab)
plt.imshow(Tab.T, cmap='gnuplot')
plt.show()
