import numpy as np

GTEMPK = 298e0
GPRES = 1e5
RH = 60
Mair = 29.06e-3
pi = np.pi
kb = 1.38064852e-23
Na = 6.02214086e+23
Rg = kb*Na

def collision_rate(diameter,mass,M):

    # viscosity of air, density oif air and mean free path in air
    viscosity     = 1.8e-5*(GTEMPK/298e0)**0.85e0  # dynamic viscosity of air
    dens_air      = Mair*GPRES/(Rg*GTEMPK)
    air_free_path = 2e0*viscosity/(GPRES*np.sqrt(8e0*Mair/(pi*Rg*GTEMPK))) # gas mean free path in air
    # knudsen number
    knudsen = 2e0 * air_free_path/diameter

    # Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34 pg 407)
    slip_correction = 1e0 + knudsen * (1.257e0 + 0.4e0*np.exp(-1.1e0/(2e0*air_free_path/diameter)))

    # particle diffusion coefficient (m^2 s^-1)
    Diff_par = (kb * GTEMPK * slip_correction) / (3e0 * pi * viscosity * diameter)

    speed_p = 0e0

    # diameter of organic compounds
    dorg = (6e0*M/Na/pi)**(1e0/3e0)               #estimated diameter (m)

    # diffusivity organic compound
    Diff_org=5e0/(16e0*Na*dorg**2e0*dens_air)*np.sqrt(Rg*GTEMPK*Mair/(2e0*pi)*((M + Mair)/M))

    speedorg=np.sqrt(8e0*kb*GTEMPK/(pi*M/Na)) #speed of organic molecules

    gasmeanfporg=3e0*(Diff_org + Diff_par)/np.sqrt(speedorg**2e0 + speed_p**2e0)

    # Knudsen number organic comp
    Knorg=2e0*gasmeanfporg/(diameter + dorg)

    # Fuchs-Sutugin correction factor for transit
    f_cororg=(0.75*(1e0 + Knorg))/(Knorg**2e0 + Knorg+0.283*Knorg + 0.75)

    Dorgeff = (Diff_org + Diff_par)*f_cororg                    # m^2/s

    collision_rate = 2e0*pi*(diameter + dorg)*Dorgeff        # mass transfer coefficient s^-1

    return collision_rate

print(collision_rate(50e-9,50e-9**3/6*np.pi*1400,np.array([200,300,4000])))
