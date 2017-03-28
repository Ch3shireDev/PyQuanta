#!/usr/bin/python
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

Nx, Ny, Nz = 2048, 256, 8
Lx, Ly, Lz = 10., 10., 1.
dt = 0.0001

dx = Lx/Nx

SpaceX = np.linspace(-Lx/2, Lx/2 , Nx)
# tab = np.exp(-SpaceX**2)
tab = np.cos(2*SpaceX/Lx*np.pi)
tab /= np.sqrt(np.linalg.norm(tab))
# tab += 100

tab = tab.astype(np.complex)
tab2 = tab.copy()
tab3 = tab.copy()

Kx = 2j*np.pi/Lx*(np.arange(Nx) - np.append(np.zeros(Nx/2), np.ones(Nx/2)*Nx))

K2 = np.exp(-1j*Kx**2*dt)

N = 5000
E = np.ones(Nx).astype(np.complex)
r = np.abs(tab2)

rr = np.abs(tab3)
ee = tab3/rr

Tab = []

E = np.ones(Nx)

for i in range(N):
    tab = ifft(K2*fft(tab))
    

    # E = np.exp(1j*np.angle(tab)) #i should remove this part somehow

    t2 = np.abs(ifft(Kx*fft(tab**2)))**2
    t3 = ifft(Kx*fft(np.abs(tab)**2))
    r4d2Theta = t2 - t3


    print i, np.amin(t2), np.amax(t2)

    # Tab.append(t2)

    plt.plot(np.abs(ifft(Kx*fft(tab)))**2) #r
    plt.show()

    #i dPsi = (dTheta**2-id2Theta)*Psi - (d2r + 2i dTheta*dR)*E
    # r = np.exp(-1j*ifft(Kx*fft(ifft(Kx*fft(E))/E))*dt)*(r-2j*ifft(Kx*fft(E))/E*ifft(Kx*fft(r))*dt)
    # E = np.exp(1j*((ifft(Kx*fft(E))/E)**2)*dt)*E
    # tab2 = r*E

    # dr = np.real(2*ifft(Kx*fft(np.abs(tab)))*ifft(Kx*fft(E))/E + np.abs(tab)*ifft(Kx**2*fft(E))/E)*dt
    # r += dr
    # RdTheta = np.real(ifft(Kx**2*fft(r)) - r*(ifft(Kx*fft(E))/E)**2)*dt
    # # tt = dr*E + 1j*E*RdTheta
    # # tab2 += E*(dr-1j*RdTheta)
    # # print np.amax(np.abs(RdTheta)), np.amax(r)



print "delta E:", np.linalg.norm(tab/np.abs(tab) - E)
print "delta R:", np.linalg.norm(np.abs(tab) - r)

plt.imshow(np.transpose(np.abs(Tab)))
# plt.imshow(np.transpose(np.angle(Tab)))
plt.show()