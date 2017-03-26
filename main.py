#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.001

dx = Lx/Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
tab = np.exp(-SpaceX**2)
tab /= np.sqrt(np.linalg.norm(tab))
tab += 10

tab = tab.astype(np.complex)
tab2 = tab.copy()

kx = 2j*np.pi/Lx*(np.arange(Nx) - np.append(np.zeros(Nx/2), np.ones(Nx/2)*Nx))

K2 = np.exp(-1j*kx**2*dt)

N = 5000
E = np.ones(Nx).astype(np.complex)
r = np.abs(tab2)

for i in range(N):
    tab = ifft(K2*fft(tab))
    E = tab/np.abs(tab) #i should remove this part somehow
    #i dPsi = (dTheta**2-id2Theta)*Psi - (d2r + 2i dTheta*dR)*E
    r = np.exp(-1j*ifft(kx*fft(ifft(kx*fft(E))/E))*dt)*(r-2j*ifft(kx*fft(E))/E*ifft(kx*fft(r))*dt)
    E = np.exp(1j*((ifft(kx*fft(E))/E)**2)*dt)*E
    tab2 = r*E

print "delta E:", np.linalg.norm(tab/np.abs(tab) - E)
print "delta R:", np.linalg.norm(np.abs(tab) - r)
print np.linalg.norm(tab - tab2)
