#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Copyright (C) 2021  Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
Contact information arca@helsinki.fi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================
"""

from pdb import set_trace as bp
# try:
import numpy as np
import sys
import os
# import json
# import requests
import pickle # optional
from scipy import optimize
# except:
    # print('To run this script, you need NumPy, SciPy, requests, pickle')

def fitSimple(T,A,B):
    return A - B/(T)

def askfor(key):
    return

def getVaps(runforever=False, args={}):
    filter_w_chem = False
    save_atoms = False
    needed = ['server','smilesfile','pram','pramfile','psat_lim','saveto']

    for i,key in enumerate(needed):
        if key in args.keys():
            pass
        else:
            return ('not enough arguments...',)
    if 'filter' in args.keys():
        if args['filter'] != []: filter_w_chem = True
    if 'save_atoms' in args.keys():
        if args['save_atoms']: save_atoms = True
    outpath, outfile = os.path.split(args['saveto'])
    ending = ''

    if args['server'] == 'AMG':
        import urllib.request
        string_a = []; count = 0
        print('Using AMG server for precalculated values')
        file = 'https://docs.google.com/spreadsheets/d/1KUfhyDydBCdQHsJF5G8pvPRZdTjcYYpsFgCPPhkj--I/export?format=csv'
        try:
            print('Contacting AMG server ... ')
            for line in urllib.request.urlopen(file):
                string = line.decode("utf-8").strip('\n')
                if not '#' in string and string.split(',')[0] != '':
                    if filter_w_chem:
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
            return('Could not contact server or got bad data, are you connected to network?',)


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
            return('Saturation vapour pressure limit was not a float.',)



    if args['server'] == 'UMan':
        try:
            f = open(args['smilesfile'])
            smilesdir, smilesfile = os.path.split(args['smilesfile'])
        except:
            print('file not found:', args['smilesfile'])
            return('file not found: '+args['smilesfile'],)


        i = 0
        Calculate_mass = False
        if args['mcm_type'] == 'csv':
            tmp_A = np.genfromtxt(args['smilesfile'], delimiter=',',dtype=None, encoding='UTF-8', unpack=True,comments='//')
            if len(tmp_A)==3:
                names, smiles, mass = tmp_A
            if len(tmp_A)==2:
                names, smiles = tmp_A
                Calculate_mass = True
                mass = np.zeros(len(names))
            smiles = np.array(smiles)
            mask = smiles != ''
            smiles = smiles[mask]
            names = np.array(names)[mask]
            mass = np.array(mass)[mask]


        if args['mcm_type'] == 'old':
            for line in f:
                if "Molecular weights for species present in the subset" in line: break
                i = i+1
            f.close()
            names, smiles, mass = np.genfromtxt(args['smilesfile'], usecols=(0,1,3),unpack=True, skip_header=i+3, dtype=str, comments='@')
        if args['mcm_type'] == 'new':
            names, smiles, mass = [],[],[]
            incomps = False
            for line in f:
                if "NameSmilesInchiInchiKeyFormulaMass" in line.replace('\t','').replace(' ',''):
                    incomps = True
                    continue
                if incomps:
                    vals = line.split()
                    if len(vals)>5:
                        names.append(vals[0])
                        smiles.append(vals[1])
                        mass.append(vals[5])
            f.close()
            names = np.array(names)
            smiles = np.array(smiles)
            mass = np.array(mass)

        N = len(names)

        n_temp = 10
        temp = np.linspace(253.15,333.15,n_temp)
        sigma = np.ones(n_temp)*0.2
        sigma[temp<268] = 0.8
        sigma[-1] = 0.8
        buffer  = np.zeros((N, n_temp))

        for i_s, cheese in enumerate(smiles):
            # get_VP(temperatures, mySmiles,vp_method='',bp_method='',umanpath='/xyz'):
            tmp = get_VP(temp,[cheese],vp_method=args['vp_method'],bp_method=args['bp_method'],umanpath=args['uMan_loc'],include_mass=Calculate_mass)
            if i_s == 0 and 'str' in str(type(tmp)):
                return ('FAILED' ,tmp)
            if i_s%100 == 0: print(f'Working with compound nr {i_s:0d}')
            buffer[i_s,:] = tmp[cheese]
            if Calculate_mass:
                mass[i_s] = tmp[cheese+'_Mw']

        n_homs = 0

        if args['pram']:
            try:
                homs = np.genfromtxt(args['pramfile'], dtype=str)
            except:
                print('PRAM file not found: "', args['pramfile']+'"')
                return('PRAM file not found: "'+args['pramfile']+'"',)

            props = os.path.split(args['pramfile'])
            props = os.path.join(props[0],props[1].replace('_names','_prop'))
            hom_mass, HomA,HomB = np.genfromtxt(props,usecols=(0,1,2), unpack=True)
            n_homs = len(homs)
            names = np.append(names,homs)
            mass = np.append(mass,hom_mass)
            try:
                hom_atoms = np.genfromtxt(props,usecols=(3,6,5,4))
                hom_atoms_full = np.zeros( (hom_atoms.shape[0],7 )  )
                hom_atoms_full[:,:4] = hom_atoms
            except:
                print('Cannot find elemental information from extra file')

        prop_matrix = np.zeros((N+n_homs, n_temp))
        prop_matrix[:N,:] = buffer

        if args['pram']:
            for i,v in enumerate(homs):
                prop_matrix[i+N,:] = fitSimple(temp,HomA[i], HomB[i] )


        a,b=np.unique(names,return_counts=True)
        if any(b>1): return ('The input produced duplicate compounds:\n'+'\n'.join(a[b>1]),)

        selected_vapours = prop_matrix[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]
        selected_vapournames = names[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]
        selected_vapourmass = mass[10**(prop_matrix[:,n_temp//2]) < float(args['psat_lim'])]
        if filter_w_chem:
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
                # print(names[j],s)
                s=s.replace('c','C')
                for i,e in enumerate(elements):
                    atoms[j,i] = (-len(s.replace(e,''))+len(s))/len(e)
                    # print(i,e,atoms[j,i])
                    s = s.replace(e,'')
                atoms[j,5] = int(round(float(mass[j]),0) - np.sum(atoms[j,:]*mass_number))

            ss      = '%-16s'+'%3i'*len(elements)
            atoms   = atoms[:,sort_ats]

            atomlib = {}
            atomlib['GENERIC'] = [22,30,0,2,0,0,0]
            for i in range(len(smiles)):
                atomlib[names[i]] = atoms[i,:]
            if args['pram']:
                for j in range(len(homs)):
                    atomlib[homs[j]] = hom_atoms_full[j,:]

        for i in range(N):
            MAB[i,0] = selected_vapourmass[i]
            try:
                MAB[i,1:], pcovS = optimize.curve_fit(fitSimple, temp, selected_vapours[i,:], maxfev = 100000, sigma=sigma)
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
        # if save_atoms: atoms = np.append(atoms, [[22,30,0,2,0,0,0]], axis=0)
    elif ii<len(names)-1 and ii>-1:
        print('GENERIC moved to last...')
        names.pop(ii)
        names.append('GENERIC')
        MAB = np.append(MAB, [MAB[ii,:]], axis=0)
        MAB = np.delete(MAB, ii, 0)
        # if save_atoms:
        #     atoms = np.append(atoms, [atoms[ii,:]], axis=0)
        #     atoms = np.delete(atoms, ii, 0)

    if save_atoms:
        path, file = os.path.split(args['saveto'])
        prefix,suffix = os.path.splitext(file)
        print('Saving elemental content to %s ...' %(os.path.join(path,prefix+'_elements'+suffix)))
        fa = open(os.path.join(path,prefix+'_elements'+suffix), 'w')
        fa.write('#Compound                    Mass               C   O   N   H   S  Cl  Br\n')

    count = 0
    f = open(args['saveto'], 'w')
    f.write( '#Compound                    Mass                        parameter_A              parameter_B\n')
    # print(names)
    for name,m,A,B in zip(names,MAB[:,0],MAB[:,1],MAB[:,2]):
        if name == 'GENERIC' or name == 'HOA':
            print('\nNOTE!: The list contains a generic/pseudovapour which is not found in the gas phase: '+name)
        if 10**(A-B/298.15)<float(args['psat_lim']):
            f.write('%-18s   %24.12f   %24.12f   %24.12f\n' %(name,m,A,B))
            count += 1
            if save_atoms:
                fa.write('%-18s   %24.12f %3d %3d %3d %3d %3d %3d %3d\n' %(name,m,*atomlib[name]))
    f.close()
    if save_atoms:
        fa.close()
    print('    Done, saved %d compounds' %count)
    return ('Done', 'Succesfully saved %d compounds' %count)


def get_VP(temperatures, mySmiles,vp_method='',bp_method='',umanpath='/xyz',include_mass=False):
    ##########################################################################################
    #											                                             #
    #    Based on an example file that loads in SMILES strings and then calculates           #
    #    pure component properties for a given temperature                                   #
    #                                                                                        #
    #    Copyright (C) 2016  David Topping : david.topping@manchester.ac.uk                  #
    #                                      : davetopp80@gmail.com     			             #
    #    Personal website: davetoppingsci.com                                                #
    #											                                             #
    #    This program is free software: you can redistribute it and/or modify                #
    #    it under the terms of the GNU Affero General Public License as published            #
    #    by the Free Software Foundation, either version 3 of the License, or                #
    #    (at your option) any later version.                                                 #
    #                                                                                        #
    #    This program is distributed in the hope that it will be useful,                     #
    #    but WITHOUT ANY WARRANTY; without even the implied warranty of                      #
    #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       #
    #    GNU Affero General Public License for more details.                                 #
    #                                                                                        #
    #    You should have received a copy of the GNU Affero General Public License            #
    #    along with this program.  If not, see <http://www.gnu.org/licenses/>.               #
    #                                                                                        #
    ##########################################################################################
    try:
        from openbabel import pybel
    except:
        print('Openbabel is needed for umansysprop. You can install it with python -m pip install openbabel')
        return('Openbabel is needed for umansysprop. You can install it with python -m pip install openbabel')
    try:
        from flask import request
    except:
        print('Flask-WTF is needed for umansysprop. You can install it with python -m pip install -U Flask-WTF')
        return('Flask-WTF is needed for umansysprop. You can install it with python -m pip install -U Flask-WTF')

    import collections

    sys.path.insert(0,os.path.abspath(umanpath))
    try:
        from umansysprop import boiling_points
        from umansysprop import vapour_pressures
        from umansysprop import critical_properties
        from umansysprop import liquid_densities
    except:
        print( 'Failed to import umansysprop, was the path to it correct?')
        return 'Failed to import umansysprop, was the path to it correct?'
    Pybel_object_dict = {}

    for Smiles in mySmiles:
       # Now create Pybel objects which are used in all property predictive techniquexs
       Pybel_object=pybel.readstring('smi',Smiles)
       Pybel_object_dict[Smiles]=Pybel_object

    ##########################################################################################
    # 2) Create a dictionary of properties based on these Pybel objects

    # NOTE: For some of the vapour pressure values, you need to perform a boiling point estimation first
    # It is therefore wise to do this initially

    # 2a) Boiling points [(K)]
    boiling_point_dict=collections.defaultdict()

    if vp_method != 'evaporation':
        for Smiles in mySmiles:
            if bp_method == 'joback_and_reid':
                boiling_point_dict[Smiles] = boiling_points.joback_and_reid(Pybel_object_dict[Smiles])
            elif bp_method == 'stein_and_brown':
                boiling_point_dict[Smiles] = boiling_points.stein_and_brown(Pybel_object_dict[Smiles])
            elif bp_method == 'nannoolal':
                boiling_point_dict[Smiles] = boiling_points.nannoolal(Pybel_object_dict[Smiles])

    vapour_pressure_dict=collections.defaultdict(list)

    for Smiles in mySmiles:
        if include_mass:
            vapour_pressure_dict[Smiles+'_Mw'] = Pybel_object_dict[Smiles].molwt
        for temperature in temperatures:
            if vp_method == 'evaporation':
                vapour_pressure_dict[Smiles].append(vapour_pressures.evaporation(Pybel_object_dict[Smiles], temperature))
            elif vp_method == 'nannoolal':
                vapour_pressure_dict[Smiles].append(vapour_pressures.nannoolal(Pybel_object_dict[Smiles], temperature, boiling_point_dict[Smiles]))
            elif  vp_method == 'myrdal_and_yalkowsky':
                vapour_pressure_dict[Smiles].append(vapour_pressures.myrdal_and_yalkowsky(Pybel_object_dict[Smiles], temperature, boiling_point_dict[Smiles]))
            else:
                return 'Unknown method '+ vp_method

    return vapour_pressure_dict



# ------------------------
if __name__ == '__main__':

    args = {}
    args['server'] = "UMan"
    args['smilesfile'] = "/home/pecl/05-ARCA/ARCA-box/ModelLib/chemistry_schemes/Full_MCM/mcm_subset_mass_reduced.txt"
    args['saveto'] = "/home/pecl/Desktop/test.txt"
    args['psat_lim'] = 1e9
    args['pram'] = False
    args['pramfile'] = ""





#
