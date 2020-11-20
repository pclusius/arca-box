import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import scipy as sc
from scipy import stats

kb  =  1.38064852e-23

create_Vap_files = False
precision = 97
suff = 'MASSIVX'
testAgainst = '../Vapour_names_98pc.dat'
limit=0
days= [4]#range(1,24)
close_figs = False

print('Aiming for %5.2f%% of total condensed mass'%precision)
hits = []
masscatch = np.zeros(len(days))
jj = -1
for day in days:
    jj +=1
    day = '%02d'%(day)
    # if len(sys.argv) > 1:
    #     file = sys.argv[1]
    # else:
    #     file = '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V100/Particles.nc'
    # try:
    #     nc = netCDF4.Dataset(Path(file), 'r')
    # except:
    #     print( 'Could not open the file, is it accessible?' )
    files = [
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V30/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V50/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V70/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V90/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V110/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V130/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V150/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V155/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V156/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V157/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V158/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V160/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V165/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V170/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V190/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V800/Particles.nc',

    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V800/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V800B100/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V800B200/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V200B40/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V200B100/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V200B200/Particles.nc',

    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/SC800/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/SC200/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/SC180/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/SC160/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/SC140/Particles.nc',

    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-10/VAERO_ALL/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-10/VAERO_70/Particles.nc',

    '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-'+day+'/OPTIFULL/Particles.nc',
    '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-'+day+'/OPTI98PC/Particles.nc',

    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/VALL50/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/VAERO_70/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V50C2/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V100C2/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V125C2/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V150C2/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V150C1/Particles.nc',
    # '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-04/V175C1/Particles.nc',
    ]

    base_set = '/home/pecl/05-ARCA/supermodel-phase-1/INOUT/HYDE/PC_2018-04-'+day+'/OPTIFULL/Particles.nc'

    for i,file in enumerate(files):
        nc = netCDF4.Dataset(Path(file), 'r')
        masses = np.sum(np.sum(nc.variables['PARTICLE_COMPOSITION'][:,:,:-2], axis=2)*nc.variables['NUMBER_CONCENTRATION'][:,:], axis=1)
        if i==0: masscatch[jj] = masses[-1]
        if i==len(files)-1: masscatch[jj] = masses[-1]/masscatch[jj]*100
        time = np.linspace(0,20,masses.shape[0])
        if i>4:
        #     plt.plot(time, masses,linestyle=':', label=labels[i])
            plt.plot(time, masses, linestyle=':', label=file[64:-12])
        else:
            plt.plot(time, masses, label=file[64:-12])
        #     plt.plot(time, masses,label=labels[i])
    # plt.ylim(0,max(masses)*1.15)
    plt.grid()
    plt.legend()
    plt.xlabel('Hours from midnight')
    plt.ylabel('kg/m3 in particles')
    if close_figs: plt.close()

    nc = netCDF4.Dataset(Path(base_set), 'r')
    cc = netCDF4.Dataset(Path(base_set[:-11]+'Chemistry.nc'), 'r')
    n_org = np.shape(nc.variables['PARTICLE_COMPOSITION'][:,:,:])[2]

    org_compo = np.sum(np.dstack([nc.variables['NUMBER_CONCENTRATION'][-1,:]]*n_org)[0,:,:]*nc.variables['PARTICLE_COMPOSITION'][-1,:,:], 0)

    plt.figure()
    org_compo[-2:] = 0
    sorter = np.argsort(-org_compo)
    compnames = netCDF4.chartostring(nc.variables['CONDENSABLES'][:])

    sortedcompo = nc.variables['PARTICLE_COMPOSITION'][:,:,:]*np.dstack([nc.variables['NUMBER_CONCENTRATION'][:,:]]*n_org)[:,:,:]
    sortedcompo[:,:,-2:] = 0
    sortedcomposum = np.sum(sortedcompo,1)[:,sorter]


    time = np.linspace(0,20,sortedcomposum.shape[0])
    maxmass = np.sum(sortedcomposum[-1,:], 0)
    first = True
    for n in range(1,len(sorter)):
        prec = (np.sum(sortedcomposum[-1,:n], 0))/maxmass*100
        if prec>=precision and first:
            n_vap_opt = n+1
            first = False
        if n%10==0:
            plt.plot(time, np.sum(sortedcomposum[:,:n], 1), label='%d first: %4.0f%%' %(n,np.sum(sortedcomposum[-1,:n], 0)/maxmass*100))
        if prec>max(99, precision):
            break
    plt.plot(time, np.sum(sortedcomposum[:,:], 1), label='All %d vapours' %(len(sorter)-2))

    plt.grid()
    plt.legend()
    plt.xlabel('Hours from midnight')
    plt.ylabel('kg/m3 in particles')
    if close_figs: plt.close()

    selecta = compnames[sorter][:n_vap_opt]
    m,a,b = np.genfromtxt('/home/pecl/05-ARCA/ChemistryPackage/Vapour_properties_all.dat', unpack=True)
    orgs = np.genfromtxt('/home/pecl/05-ARCA/ChemistryPackage/Vapour_names_all.dat', dtype=str)

    GENindex = list(orgs).index('GENERIC')

    if create_Vap_files:
        opt_names = open('../Vapour_names_%s_%s.dat'%(suff,day),'w')
        opt_props = open('../Vapour_properties_%s_%s.dat'%(suff,day),'w')

    plt.figure()
    plt.title('%d vapours'%n_vap_opt)

    names = []
    testAgainstNames = np.genfromtxt(testAgainst, dtype=str, skip_footer=1)
    if limit>0:testAgainstNames = testAgainstNames[:limit]
    for i in range(n_vap_opt):
        names.append(selecta[i].strip())
        if create_Vap_files: opt_names.write(selecta[i]+'\n')
        conc = np.mean(cc.variables[selecta[i].strip()][50:])
        ind = list(orgs).index(selecta[i].strip())
        plt.loglog(sortedcomposum[-1,i], conc/10**(a[ind]-b[ind]/298), 'ro')
        plt.text(sortedcomposum[-1,i], conc/10**(a[ind]-b[ind]/298),selecta[i].strip(),  fontsize=6)
        if create_Vap_files: opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[ind],a[ind],b[ind]))

    if create_Vap_files:
        opt_names.write('GENERIC\n')
        opt_props.write('%24.12f\t%24.12f\t%24.12f\n' %(m[GENindex],a[GENindex],b[GENindex]))

        opt_names.close()
        opt_props.close()

    common = len(list(set(names).intersection(testAgainstNames)))
    nt = len(testAgainstNames)
    na = len(names)
    hits.append(common/na*100)
    print(day+':',common,'vapours of', nt,round(common/nt*100,0), '% exist in actual (total', na,round(common/na*100,0),'%)', 'Total mass captured:%4.0f%%'%masscatch[jj])
    plt.xlabel('pp mass in the end')
    plt.ylabel(r'$\overline{C}$/P*')
    if close_figs: plt.close()

plt.figure()
plt.bar(days,masscatch)
plt.title('Goodness of current vap file: '+testAgainst)
plt.xlabel('day')
plt.ylabel('% of total mass modelled')
print('On average %5.1f%%'%np.mean(hits), 'of correct vapours catched')
print('On average %5.1f%%'%np.mean(masscatch), 'of total mass catched')
plt.show()



#
