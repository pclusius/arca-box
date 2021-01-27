import numpy as np
import re
import sys
import os

def multi(x, i=-1):
    return x

def shift(x, i=-1):
    return x

batch = None
mod = -1

if len(sys.argv)==1:
    print('Please give your batch file path as input (e.g. ../INOUT/PC_CONTROL_0000-0004.bash)')

else:
    batch = sys.argv[1]
    batchdir = os.path.split(batch)[0]
    mod = int(sys.argv[2])

if batch != None:

    if os.path.exists(batch):
        initfile = []
        with open(batch) as f:
            for line in f:
                vars = line.split()
                if 'cd' in vars[0]: cd = vars[1]
                if 'arcabox.exe' in vars[0]:
                    initfile.append (os.path.join(batchdir,cd,vars[1]))


        for i,init in enumerate(initfile):
            found = False
            if os.path.exists(init): os.rename(init, init+'.bkup')
            fw = open(init, 'w')
            with open(init+'.bkup') as f:
                for line in f:
                    vars = line.split()
                    if len(vars) > 0:
                        n = r'MODS\(0*%d\)'%mod
                        if re.findall(n,vars[0]):
                            found = True
                            vars[4] = ('%12.8e'%(multi(float(vars[4].replace('d', 'e')), i))).replace('e','d')
                            vars[5] = ('%12.8e'%(shift(float(vars[5].replace('d', 'e')), i))).replace('e','d')
                            print('Found index %d from file '%mod+init)
                    fw.write(' '.join(vars)+'\n')
                if not found: print('index %d was not found from file '%mod+init)

    else:
        print('Could not find the file')
