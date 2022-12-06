import numpy as np
import re
import sys
import os
import shutil
from pdb import set_trace as bp

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

# ---------------------------------------------------------------------------------------------------------------------
# ----------- USER OPTIONS START HERE ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

"""
This matrix defines the variations. Number of steps is used for both multi and shift, usually one or the other stays
constant, meaning that both min and max are the same. Steps are taken linearly if "log" is not 1. Call the script with
the path to the bash-script that was created with ARCA's batch tool, e.g.:

python3 Variations.py ../../INOUT/PC_CONTROL_2018-04-01-2018-04-05.bash

The script will first test that the indices are found in all relevant files, the you are asked if you want to continue
with file operations and modifying the bash file. if you type y, N new initfiles and empty target directories will be
created, where N = product of values in column "number of steps".
"""
ops = np.array([
# index    number-of-steps   min_multi   max_multi    min_shift    max_shift     log
# [ 1 ,             2 ,            1 ,         1 ,          0  ,        0   ,       0],
# [ 2 ,             4 ,            1 ,         1 ,          0  ,        0   ,       0],
# [ 3 ,             3 ,            1 ,         1 ,          0  ,        0   ,       0],

[ [1,2,3]  ,       5  ,            0.8 ,      1.2 ,         0 ,         0 ,         0 ],
[ 4        ,       5  ,            0.8 ,      1.2 ,         0 ,         0 ,         0 ],
],dtype=object)

# ---------------------------------------------------------------------------------------------------------------------
# ----------- USER OPTIONS END HERE -----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

def readfile(init,multi,shift,mod,dir=None):
    with open(init) as f:
        modified_file = []
        found = False
        rundirs = ''
        indices = {}
        group=False
        if 'list' in str(type(mod)):
            n = [r'MODS\(0*%d\)'%m for m in mod]
            group=True
        else:
            n = [r'MODS\(0*%d\)'%mod]
        c = r'RUN_NAME'
        # bp()
        for line in f:
            vars = line.split()
            if len(vars) > 0:
                if re.findall(c,vars[0]):
                    rundirs = vars[2].replace("'", "")
                    if dir != None: vars[2] = "'"+dir+"'"
                if re.findall(r'MODS\((\d+)\)', vars[0]):
                    iii = (int(re.findall(r'MODS\((\d+)\)', vars[0])[0]))
                    indices[vars[-1]] = iii
                if any (re.findall(nn,vars[0]) for nn in n):
                    # bp()
                    found = True
                    vars[4] = ('%12.8e'%(multi*float(vars[4].replace('d','e')))).replace('e','d')
                    vars[5] = ('%12.8e'%(shift+float(vars[5].replace('d','e')))).replace('e','d')
                    print('Found %s from file '%vars[0]+init)
                if vars[0][0] != '&' and vars[0][0] != '#' and vars[0][0] != '/':
                    vars[0] = ' '+vars[0]
            modified_file.append(' '.join(vars)+'\n')

    return modified_file, rundirs, found, indices


def zzzz(batch, batchdir, ops, dryrun=False, nopause=False):
    if batch == None:
        return 'No input'
    if not os.path.exists(batch):
        return 'Could not find the file'

    initfile = []
    batchfile = []
    bashhdr = []
    with open(batch) as f:
        for line in f:
            vars = line.split()
            if len(vars)>0:
                if 'cd' in vars[0]: cd = vars[1]
                if 'arcabox.exe' in vars[0]:
                    initfile.append (os.path.join(batchdir,cd,vars[1]))
                    batchfile.append (vars)
                else:
                    bashhdr.append(line.strip('\n'))
    rundirs = []
    newfiles = []
    dim = ops.shape[0]
    N = int(np.prod(ops[:,1]))
    e = int(np.log10(N))+2
    mout = np.zeros((N,dim))
    sout = np.zeros((N,dim))
    M = []
    S = []
    # bp()
    for n in range(dim):
        if ops[n,6] == 1:
            M.append(np.logspace(np.log10(ops[n,2]),np.log10(ops[n,3]),int(ops[n,1])))
            S.append(np.logspace(np.log10(ops[n,4]),np.log10(ops[n,5]),int(ops[n,1])))
        else:
            M.append(np.linspace(ops[n,2],ops[n,3],int(ops[n,1])))
            S.append(np.linspace(ops[n,4],ops[n,5],int(ops[n,1])))
    # bp()
    for i in range(N):
        d = 1
        for j in range(dim):
                # bp()
                mout[i,j] = M[j][(i//d)%len(M[j])]
                sout[i,j] = S[j][(i//d)%len(S[j])]
                d = d*len(M[j])


    for r in ['test', 'forreal']:
        for i_init,init in enumerate(initfile):
            if not os.path.exists(init): return 'Initifile %s not found' %init
            if r == 'forreal':
                print('\nNow working with ',init,'\n')

                # Here the "old" bash file row is appended to the new bash file text
                newfiles.append(' '.join(batchfile[i_init]))

                # Loop through all the initfiles
                for fn in range(N):

                    # formatting for the index number
                    format = '%%0%dd_'%e

                    # Pick the initfile name from the bash command, append index number
                    newfile = os.path.join(os.path.split(init)[0], format%(1+fn)+os.path.split(init)[1])

                    # Create the new initfile, at this point it is just a copy of the original
                    shutil.copyfile(init, newfile)

                    # workdir is the RUN_NAME from initfile
                    workdir = format%(1+fn)+rundirs[i_init]

                    # parse the new line for bash script
                    to_bashfile = (
                        batchfile[i_init][0]  # the ./arcabox.exe part
                        +' '+ os.path.join(os.path.split(batchfile[i_init][1])[0], os.path.split(newfile)[1]) # new initfile
                        +' '+ batchfile[i_init][2] # the tee command
                        +' '+ os.path.join(os.path.split(os.path.split(batchfile[i_init][3])[0])[0],workdir,'runReport.txt')  # the runReport path
                        )
                    # save it for later...
                    newfiles.append(to_bashfile)

                    # If the output directories did not exist, create them
                    if not os.path.exists(os.path.join(os.path.split(init)[0],'..',workdir)):
                        os.mkdir(os.path.join(os.path.split(init)[0],'..',workdir))

                    # Now the silly part, open and write the files N times. Inefficient but won't be the bottleneck in your work...
                    for dd in range(dim):
                        a,_,_,_ = readfile(newfile, mout[fn,dd],sout[fn,dd],ops[dd,0],workdir)
                        fff = open(newfile,'w')
                        fff.write(''.join(a))
                        fff.close()

            else:
                for i_c,mod in enumerate(ops[:,0]):
                    # bp()
                    modified_file, rd, found, indices = readfile(init,1.0,0.0,mod)
                    if dryrun:
                        return indices
                    rundirs.append(rd)
                    if not found:
                        return 'index %d was not found from file '%mod+init

        if r == 'forreal':
            bash = bashhdr
            bash += ['','#  ind  steps mul_mn          mul_mx          sh_mn           sh_mx           log']
            for i in range(ops.shape[0]):
                if 'list' in str((ops[i,:])):
                    buff = [','.join(['%d'%indzz for indzz in ops[i,0]])]
                    # bp()
                    buff += list(ops[i,1:])
                    bash += ['#  '+'{:s}    {:0.0f}     {:12.8e}  {:12.8e}  {:12.8e}  {:12.8e}  {:0.0f}'.format(*buff) ] + ['']
                else:
                    buff = ops[i,:]
                    bash += ['#  '+'{:0.0f}    {:0.0f}     {:12.8e}  {:12.8e}  {:12.8e}  {:12.8e}  {:0.0f}'.format(*buff) ] + ['']
            bash += newfiles

            bf = batch.replace('.bash','')+'_mod'+'.bash'
            f = open(bf, 'w')
            f.write('\n'.join(bash))
            f.close()
            return '\nCompleted. The new bash file is %s\n\nThank you for flying with us!\n'%bf

        elif not nopause:
            cont = input('continue? \n')
            if cont.upper() == 'Y' or cont.upper() == 'YES':
                pass
            else: return 'Okay then. Bye.'

if __name__ == '__main__':
    if len(sys.argv)==1:
        print('Please give your batch file path as input (e.g. ../INOUT/PC_CONTROL_0000-0004.bash)')

    else:
        nopause = False
        dryrun  = False
        batch = sys.argv[1]
        if len(sys.argv)>2:
            if sys.argv[2] == 'nopause':
                nopause = True
            if sys.argv[2] == 'dryrun':
                dryrun = True
        batchdir = os.path.split(batch)[0]

        print(zzzz(batch, batchdir, ops))
