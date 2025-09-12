import numpy as np
import os,sys,netCDF4
import matplotlib.pyplot as plt

spectra = False
writeFiles = False
readNames = False
pathToBinaryFile = ''
selectTime = 9999999999
chemistryPath = ''


def read_bin(file='TIMESERIES.r16'):
    f = open(file, 'rb')
    conc = np.fromfile(f, dtype=np.float64)
    dt = conc[0]
    d0 = int(conc[1])
    d1 = int(conc[2])
    conc = conc[3:].reshape([d0,d1], order='F')
    return dt,d0,d1,conc

def get_conc_at_time(file, time=None):
    dt,d0,d1, conc = read_bin(file)
    full_time = (d1-1)*dt
    out = np.zeros((d0,2))
    # pick the nerest index
    if time is None or float(time)>=full_time:
        i_time = d1-1
    else:
        i_time = int(round(float(time)/dt))
    out[:,0] = conc[:,i_time]
    if i_time+1 == conc.shape[1]:
        out[:,1] = (conc[:,i_time] - conc[:,i_time-1])/dt
    else:
        out[:,1] = (conc[:,i_time+1] - conc[:,i_time])/dt
    return dt,d0,d1,out

def get_names(file,chemistryPath):
    nc = netCDF4.Dataset(file+'/Chemistry.nc')
    names = []
    inds = []
    chem = nc.Chemistry_module
    if chemistryPath != '':
        path = os.path.join(chemistryPath,'Elements.txt')
    else:
        path = os.path.join('src','chemistry',chem,'Elements.txt')
    comnames = np.genfromtxt(path, usecols=0, dtype=str)
    mass = np.genfromtxt(path, usecols=1)
    C_H_O_N = np.genfromtxt(path, usecols=(2,3,4,5),dtype=int)
    rad = np.genfromtxt(path, usecols=10,dtype=int)
    cn = {comnames[i]:C_H_O_N[i,:] for i in range(len(comnames))}
    rd = {comnames[i]:rad[i] for i in range(len(comnames))}
    ms = {comnames[i]:mass[i] for i in range(len(comnames))}
    for k in nc.variables.keys():
        try:
            if nc.variables[k].type == 'gas':
                names.append(k)
                inds.append(nc.variables[k].index-1)
        except:
            pass
    return list(np.array(names)[sorted(inds)]), cn, rd, ms


def update(selectTime=selectTime,
    pathToBinaryFile=pathToBinaryFile,
    readNames=readNames,
    writeFiles=writeFiles,chemistryPath=chemistryPath):
    if readNames:
        names,cn,rd, mass = get_names(os.path.split(pathToBinaryFile)[0],chemistryPath)
        M = np.array([mass[n] for n in names])
    dt,nComp,nTime,Conc = get_conc_at_time(pathToBinaryFile,selectTime)
    if readNames:
        RO2_pool = np.sum([Conc[i_n,0]*rd[n] for i_n,n in enumerate(names)])
        NO, HO2 = [Conc[names.index(n),0] for n in ['NO','HO2']]
    if selectTime>dt*nTime: selectTime = dt*(nTime-1)
    if writeFiles:
        np.savetxt(os.path.split(pathToBinaryFile)[0]+'/CdC.txt', Conc, header=f'Concentrations and change of concentrations at time {selectTime} s, molec cm-3 (s-1)')
        if readNames:
            with open(os.path.split(pathToBinaryFile)[0]+'/Names.txt', 'w') as file:
                file.write('#-name--------------------------C--H--O--N-\n')
                for i_n, name in enumerate(names):
                    atoms = cn.get(name, [0,0,0,0])
                    file.write(f'{name:30s} {atoms[0]:2d} {atoms[1]:2d} {atoms[2]:2d} {atoms[3]:2d}\n')

            with open(os.path.split(pathToBinaryFile)[0]+'/NO_HO2_sumRO2.txt', 'w') as file:
                file.write(f'{NO} {HO2} {RO2_pool}\n')
    return Conc, names, M, RO2_pool,NO,HO2

def updatePlot(C,M,ax,mlim=100,w=0.4):
    ax.bar(M[M>mlim], C[M>mlim,0], width=w)
    ax.set_yscale('log')
    ax.set_xlabel('Mass (average molar, Da)')
    ax.set_ylabel('cm$^{-3}$')
    # plt.title(f'Concentrations at {selectTime:0.2f} s')

if __name__ == '__main__':
    for ia, arg in enumerate(sys.argv):
        if arg=='-f':
            writeFiles = True
        if arg=='-n':
            readNames = True
        if arg=='-t':
            selectTime = float(sys.argv[ia+1])
        if arg=='-p':
            pathToBinaryFile = sys.argv[ia+1]
        if arg=='-c':
            chemistryPath = sys.argv[ia+1]
        if arg=='-s':
            spectra = True
    if pathToBinaryFile == '':
        exit()
    print(f"""
Selected options:
writeFiles (-f) = {writeFiles}
readNames (-n) = {readNames}
selectTime (-t) = {selectTime}
pathToBinaryFile (-p) = {pathToBinaryFile} the binary file (e.g. TIMESERIES.r8)
chemistryPath (-c) = {chemistryPath} set path to chemistry if this tool is not called from ARCA root
spectra (-s) = {spectra} if set, will plot mass spectra
    """)

    # print("To read again data, edit selectTime and run:")
    # print("Conc, names,M, RO2_pool,NO,HO2 = update(selectTime=selectTime,pathToBinaryFile=pathToBinaryFile,readNames=readNames,writeFiles=writeFiles,chemistryPath=chemistryPath)")

    Conc, names,M, RO2_pool,NO,HO2 = update(selectTime=selectTime,
        pathToBinaryFile=pathToBinaryFile,
        readNames=readNames,
        writeFiles=writeFiles,chemistryPath=chemistryPath)
    # plt.ion()
    if spectra:
        f,ax=plt.subplots(1, figsize=(8,4))
        # print('To update plot, run: updatePlot(Conc,M,ax=ax)')
        updatePlot(Conc,M,ax=ax)
        plt.show()
#
