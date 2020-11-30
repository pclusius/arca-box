import netCDF4
import numpy as np
import sys

GTEMPK = 298e0
GPRES = 1e5
RH = 60
Mair = 29.06e-3
pi = np.pi
kb = 1.38064852e-23
Na = 6.02214086e+23
Rg = kb*Na

if len(sys.argv) > 1:
    N_vapours = sys.argv[1]
else:
    N_vapours=1861

if len(sys.argv) > 2:
    file = sys.argv[2]
else:
    file = '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/VALL50/Chemistry.nc'

opt_names = open('../Vapour_names_C1.dat','w')
opt_props = open('../Vapour_properties_C1.dat','w')


def criteria(aa, bb, M, method):
    if method==1:
        return aa/bb
    if method==2:
        return aa/bb*(collision_rate(50e-9,50e-9**3/6*np.pi*1400,M))
    if method==3:
        f = aa/bb
        f = np.where(f>0, np.log(f), -np.inf())
        print(f)
        return f

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
    speed_p = np.sqrt(8E0*kb*GTEMPK/(pi*mass))
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
# print(collision_rate(50e-9,50e-9**3/6*np.pi*1400,np.array([200,300,4000])))


try:
    nc = netCDF4.Dataset(file, 'r')
except:
    print( 'Could not open the file, is it accessible?' )

C_all = {}
for comp in nc.variables.keys():
    C_all[comp] = np.mean(nc.variables[comp][50:])

m,a,b = np.genfromtxt('/home/pecl/05-ARCA/ChemistryPackage/Vapour_properties_all.dat', unpack=True)
orgs = np.genfromtxt('/home/pecl/05-ARCA/ChemistryPackage/Vapour_names_all.dat', dtype=str)
conc_sat = 10**(a-b/298.)

GENindex = list(orgs).index('GENERIC')
testvalues = np.zeros(len(orgs))

for i,comp in enumerate(orgs):
    if i != GENindex:
        # testvalues[i] = criteria(C_all[comp], conc_sat[i],m[i]*1e-3,3)
        ff = C_all[comp]**2/conc_sat[i] # C2
        ff = C_all[comp]/conc_sat[i] # C1
        if ff>0:
            f = np.log(ff)
        else:
            f = -np.inf
        print(comp, f, conc_sat[i])
        testvalues[i] = f



order = np.argsort(-testvalues)

for i in range(N_vapours):
    # print(orgs[order][i], testvalues[order][i])
    if 'CH4' == orgs[order][i]:
        print('##############################################3')
    opt_names.write(orgs[order][i]+'\n')
    opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[order][i],a[order][i],b[order][i]))

opt_names.write('GENERIC\n')
opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[GENindex],a[GENindex],b[GENindex]))

opt_names.close()
opt_props.close()




#
