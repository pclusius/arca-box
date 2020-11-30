#MassDiffusivity.py

import numpy as np
import matplotlib.pyplot as plt
Rg   =  8.31446
Na  =  6.02214086e+23
Omega = 3.1

def DChapman(MX, VolX, Mair=28.98, p=1.01325e5, T=293):
    return 1e-7*T**(1.75)*np.sqrt(1/(Mair)+1/(MX))/(p/1.01325e5*0.25*((VolX*6/np.pi)**(1e0/3e0)+20.1**(1e0/3e0))**2e0 * Omega)

def DARCA(MX, VolX, Mair=28.98, p=1.01325e5, T=293):
    return 1e-7*T**(1.75e0)*np.sqrt(1/(Mair)+1/(MX))/(p/1.01325e5*(VolX**(1e0/3e0)+20.1**(1e0/3e0))**2e0)

def DIFF(MX, dia, T=293, p=1.01325e5, Mair=28.98):
    Mair = Mair*1e-3
    MX = MX*1e-3
    return 5e0/(16e0*Na* dia**2e0*(Mair * p/(Rg*T))) * np.sqrt(Rg*T*Mair/(2.*np.pi)*((MX+Mair)/MX))


def vol(NrC,NrO,NrH,NrN,NrS):
    return NrC*15.9e0+NrO*6.11e0+NrH*2.31e0+NrN*4.54e0+NrS*22.9e0

volSA = vol(0,4,2,0,1)
volBZ = vol(6,0,0,6,0)
volC2H5OH = vol(2,1,0,6,0)
propanol = vol(3,1,0,8,0)

vol=volSA


print('DChapman:',DChapman(98, vol))
print('DARCA   :',DARCA(98, vol))

T = np.linspace(260,330,100)
p = np.linspace(700,1050,100)*1e2

plt.figure()
plt.plot(T, DChapman(98, vol, T=T), label='Chapman')
plt.plot(T, DARCA(98, vol, T=T), label='ARCA')
plt.legend()

# plt.figure()
# plt.plot(p, DChapman(98, vol, p=p), label='Chapman')
# plt.plot(p, DARCA(98, vol, p=p), label='ARCA')
# plt.legend()

plt.show()
