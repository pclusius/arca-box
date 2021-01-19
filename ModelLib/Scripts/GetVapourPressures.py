#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Petri CLusius
# Created Date: 18/01/2021
# =============================================================================
try:
    import numpy as np
    import sys
    import os
    import json
    import requests
    import pickle # optional
    from scipy import optimize
except:
    print('To run this script, you need NumPy, SciPy, requests, pickle')

def fitSimple(T,A,B):
    return A - B/(T)

def ManU(temp, smiles):
    N=len(smiles)
    prop_matrix = np.zeros((N, len(temp) ) )
    url = 'http://umansysprop.seaes.manchester.ac.uk/api/vapour_pressure'
    print('Contacting ', url)
    header = {
            'Accept': 'application/json',
            'Content-Type': 'application/json',
            'Accept-Encoding': 'gzip, deflate',
            'User-Agent': 'python-requests/2.4.1 CPython/2.7.6 Linux/3.13.0-48-generic',
            'keep-alive': 'timeout=600, max=100',
            'Connection': 'keep-alive'
            }

    get_vapours = {
                "vp_method": "nannoolal",
                "bp_method": "nannoolal",
                "temperatures": list(temp),
                "compounds": list(smiles)
                }

    x = requests.post(url, data = json.dumps(get_vapours), headers=header)
    data = x.json()[0]
    alldata = data['data']
    # print(len(alldata))
    for i in range(len(alldata)):
        prop_matrix[i%N,i//N] = float(alldata[i]['value'])
    return prop_matrix

def askfor(key):
    return

def getVaps(runforever=False, slave=True, args={}):
    filter = False
    save_atoms = False
    needed = ['server','smilesfile','pram','pramfile','psat_lim','saveto']

    for i,key in enumerate(needed):
        if key in args.keys():
            pass
        else:
            if not slave:
                args[key] = askfor(key)
            else:
                print('not enough arguments...')
                return
    if 'filter' in args.keys():
        if args['filter'] != []: filter = True
    if 'save_atoms' in args.keys():
        if args['save_atoms']: save_atoms = True
    outpath, outfile = os.path.split(args['saveto'])
    ending = ''

    if args['server'] == 'AMG':
        import urllib.request
        string_a = []; count = 0
        print('Using AMG server for precalculated values')
        file = 'https://docs.google.com/spreadsheets/d/1i1vs2LeJkj-fwfo-294fejxGipcJ1IPCQi0bIQT5MPY/export?format=csv'
        try:
            print('Contacting AMG server ... ')
            for line in urllib.request.urlopen(file):
                string = line.decode("utf-8").strip('\n')
                if not '#' in string and string.split(',')[0] != '':
                    if filter:
                        if string.split(',')[0].upper().strip() in args['filter']:
                            string_a.append(string.split(','))
                            count += 1
                            ending = ', filtered with current chemistry.'
                    else:
                        string_a.append(string.split(','))
                        count += 1
            print('    Retrieved %d compounds %s'%(count,ending))
        except:
            print('Could not contact server or got bad data, are you connected to network?')
            return

        N = len(string_a)
        MAB = np.zeros((N,3))
        atoms = np.zeros((N,7))
        names = [(string_a[i][0]) for i in range(N)]
        MAB[:,0] = [float(string_a[i][1]) for i in range(N)]
        MAB[:,1] = [float(string_a[i][2]) for i in range(N)]
        MAB[:,2] = [float(string_a[i][3]) for i in range(N)]
        for j in range(7):
            atoms[:,j] = [float(string_a[i][5+j]) for i in range(N)]
        try:
            print('Saving compounds with Psat < %6.2e to %s ... '%(float(args['psat_lim']), args['saveto']))
        except:
            print('Saturation vapour pressure limit was not a float.')
            return


    if args['server'] == 'UMan':
        try:
            f = open(args['smilesfile'])
            smilesdir, smilesfile = os.path.split(args['smilesfile'])
        except:
            print('file not found:', args['smilesfile'])
            return

        i = 0
        for line in f:
            if "Molecular weights for species present in the subset" in line: break
            i = i+1
        f.close()

        names, smiles, mass = np.genfromtxt(args['smilesfile'], usecols=(0,1,3),unpack=True, skip_header=i+3, dtype=str, comments='@')
        N = len(names)

        n_temp = 8
        temp = np.linspace(253.15,323.15,n_temp)

        n_first = (len(names)%100)
        if n_first<len(names):
            rest = len(names[n_first:])//100

        buffer  = np.zeros((N, n_temp))
        try:
            buffer = pickle.load(open(os.path.join(smilesdir, smilesfile+"_UMan_Fetch.pickle"), "rb"))
        except:
            if rest>0:
                print('Fetching compounds in chunks of 100, have patience, this might take a while.')
                for j in range(rest):
                    print('Fetching compounds from %i to %i of %d.'%(j*100, (j+1)*100, N))
                    buffer[j*100:(j+1)*100,:] = ManU(temp, smiles[j*100:(j+1)*100])
            print('Fetching compounds %i to %i.'%(N-n_first+1, N))
            buffer[-n_first:,:] = ManU(temp, smiles[-n_first:])
            pickle.dump( buffer, open( os.path.join(smilesdir, smilesfile+"_UMan_Fetch.pickle"), "wb" ) )

        n_homs = 0

        if args['pram']:
            homs = np.genfromtxt(args['pramfile'], dtype=str)
            props = os.path.split(args['pramfile'])
            props = os.path.join(props[0],props[1].replace('_names','_prop'))
            hom_mass, HomA,HomB = np.genfromtxt(props,usecols=(0,1,2), unpack=True)
            n_homs = len(homs)
            names = np.append(names,homs)
            mass = np.append(mass,hom_mass)
            try:
                hom_atoms = np.genfromtxt(props,usecols=(3,4,5,6))
            except:
                print('Cannot find elemental information from extra file')

        prop_matrix = np.zeros((N+n_homs, n_temp))
        prop_matrix[:N,:] = buffer

        if args['pram']:
            for i,v in enumerate(homs):
                prop_matrix[i+N,:] = fitSimple(temp,HomA[i], HomB[i] )

        selected_vapours = prop_matrix[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]
        selected_vapournames = names[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]
        selected_vapourmass = mass[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]

        if filter:
            filter_a = [False]*len(selected_vapours)
            for i in range(len(selected_vapours)):
                if selected_vapournames[i].upper() in args['filter']:
                    filter_a[i] = True
                    ending = ', filtered with current chemistry.'

            selected_vapours = selected_vapours[filter_a]
            selected_vapournames = selected_vapournames[filter_a]
            selected_vapourmass = selected_vapourmass[filter_a]

        N = np.size(selected_vapours,0)
        print(N, ' compounds selected'+ending)
        MAB = np.zeros((N,3))
        # np.savetxt('Vapour_names.dat', selected_vapournames, fmt='%s')

        # Determine elemental composition from SMILES
        elements    = np.array(['Br','Cl','C','O','N','H','S'], dtype=str)
        atoms       = np.zeros((len(smiles)+n_homs,len(elements)), dtype=int)
        if save_atoms:
            #                        0    1    2   3   4   5   6
            mass_number = np.array([80,35,12,16,14,1,32])
            sort_ats    = [2,3,4,5,6,1,0]

            for j,s in enumerate(smiles):
                for i,e in enumerate(elements):
                    atoms[j,i] = (-len(s.replace(e,''))+len(s))/len(e)
                    s = s.replace(e,'')
                atoms[j,5] = int(round(float(mass[j]),0) - np.sum(atoms[j,:]*mass_number))

            ss      = '%-16s'+'%3i'*len(elements)
            atoms   = atoms[:,sort_ats]

            atomlib = {}
            for i in range(len(smiles)):
                atomlib[names[i]] = atoms[i,:]
            if args['pram']:
                for j in range(len(homs)):
                    atomlib[homs[j]] = hom_atoms[j,:]

        for i in range(N):
            MAB[i,0] = selected_vapourmass[i]
            try:
                MAB[i,1:], pcovS = optimize.curve_fit(fitSimple, temp, selected_vapours[i,:], maxfev = 100000)
            except:
                print('WARNING: Failed to fit Antoine equation to : '+selected_vapournames[i])

        names = selected_vapournames.tolist()

    ii = -1
    if 'HOA' in names:
        names[names.index('HOA')] = 'GENERIC'
        print('Changed "HOA" to "GENERIC" as is the current convention')
    if 'GENERIC' in names:
        ii = names.index('GENERIC')
    if ii == -1:
        names.append('GENERIC')
        MAB = np.append(MAB, [[437,10,10000]], axis=0)
        if save_atoms: atoms = np.append(atoms, [[22,30,0,2,0,0,0]], axis=0)
    elif ii<len(names)-1 and ii>-1:
        print('GENERIC moved to last...')
        names.pop(ii)
        names.append('GENERIC')
        MAB = np.append(MAB, [MAB[ii,:]], axis=0)
        MAB = np.delete(MAB, ii, 0)
        if save_atoms:
            atoms = np.append(atoms, [atoms[ii,:]], axis=0)
            atoms = np.delete(atoms, ii, 0)

    if save_atoms:
        path, file = os.path.split(args['saveto'])
        print('Saving elemental content to %s ...' %(os.path.join(path,'elements_'+file)))
        fa = open(os.path.join(path,'elements_'+file), 'w')
        fa.write('#Compound                    Mass               C   O   N   H   S  Cl  Br\n')

    count = 0
    f = open(args['saveto'], 'w')
    f.write( '#Compound                    Mass                        parameter_A              parameter_B\n')
    for name,m,A,B,ats in zip(names,MAB[:,0],MAB[:,1],MAB[:,2],atoms):
        if name == 'GENERIC' or name == 'HOA':
            print('\nNOTE!: The list contains a generic/pseudovapour which is not found in the gas phase: '+name)
        if 10**(A-B/298.15)<float(args['psat_lim']):
            f.write('%-18s   %24.12f   %24.12f   %24.12f\n' %(name,m,A,B))
            count += 1
            if save_atoms:
                fa.write('%-18s   %24.12f %3d %3d %3d %3d %3d %3d %3d\n' %(name,m,*ats))
    f.close()
    if save_atoms:
        fa.close()
    print('    Done, saved %d compounds' %count)


# ------------------------
if __name__ == '__main__':
    # if len(sys.args)
    root = '/home/pecl/05-ARCA/ChemistryPackage/mcm_large.txt'
    qserver = input('Use AMG database y/n (default=yes)?:\n')
    if qserver.upper()=='N':
        server = False
        mcmfile = input('Give path to SMILES file (from (MCM), default %s):\n'%root)
        if mcmfile != '': root = mcmfile

        replace_UMAN = input('Fetch from UMAN database y/n (replacing any current, default=No)?:\n')
        if replace_UMAN.upper() == 'Y' and os.path.exists('umansysprop.pickle'):
            os.remove('umansysprop.pickle')
    else: server = True

    qpram = input('Use PRAM y/n (default=yes)?:\n')
    if qpram.upper()=='N':
        pram = False
    else: pram = True

    psat_lim = input('Limit Psat[atm] to (default=1e-6)?:\n')
    if psat_lim=='':
        psat_lim = 1e-6
    else:
        psat_lim = float(psat_lim)

    qsortDesc = input('Sort to descending order by Psat (y/n, default = No sorting)?:\n')
    if qsortDesc.upper()=='Y':
        sortDesc = True
    else: sortDesc = False

    print('Using ',root)
    print('Picking compounds with vapour pressure lower than ',psat_lim)


    runon = True
    while runon:
        runon = False
        getVaps(runforever=True, slave=False)
        q = input('\nRestart (r) program or quit (q, default)?:\n')
        if q.lower() == 'r': runon = True

#
