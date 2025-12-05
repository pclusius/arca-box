import numpy as np

def collision_rate(diameter, mass,temp,press):
    if 'NoneType' in str(type(mass)):
        mass = 1400*np.pi*diameter**3/6
    Rg = 8.31446
    kb = 1.38064852e-23
    Mair = 28.966000e-3
    M = 98.07e-3
    thermal_speed = np.sqrt(8e0*kb*temp/(np.pi*M/(Rg/kb)))

    diff_vol = 4*6.11 + 2*2.31 + 22.9
    diffusivity = 1e-7 * temp**1.75 * np.sqrt( 1/(Mair*1e3) + 1/(1e3*M) ) / (press/1.01325e5 * (diff_vol**(1./3.) + 20.1**(1./3.))**2)
    diffusivity_dia = (6 * (1e-30 * diff_vol) / np.pi )**(1/3.)
    # viscosity of air and mean free path in air
    viscosity     = 1.8e-5*(temp/298e0)**0.85e0  # dynamic viscosity of air
    air_free_path = 2e0*viscosity/(press*np.sqrt(8e0*Mair/(np.pi*Rg*temp))) # gas mean free path in air

    # Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34 pg 407)
    slip_correction = 1e0 + 2e0*air_free_path/diameter * (1.257e0 + 0.4e0*np.exp(-1.1e0/(2e0*air_free_path/diameter)))

    Diff_par = (kb * temp * slip_correction) / (3e0 * np.pi * viscosity * diameter)

    speed_p = np.sqrt(8e0*kb*temp/(np.pi*mass))

    vapmeanfp = 3e0*(diffusivity + Diff_par)/np.sqrt(thermal_speed**2 + speed_p**2)

    Knudsen = 2e0 * vapmeanfp/(diameter + diffusivity_dia)

    # Fuchs-Sutugin correction factor for transition regime
    FS_corr        = (0.75 * (1e0 + Knudsen)) /(Knudsen**2 + Knudsen + 0.283 * Knudsen + 0.75e0)
    D_vap_eff      = (diffusivity + Diff_par) * FS_corr # m^2/s
    collision_rate = 2e0*np.pi*(diameter + diffusivity_dia) * D_vap_eff       # mass transfer coefficient s^-1
    return collision_rate
