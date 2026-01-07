import numpy as np
import os, re, pickle, gzip,netCDF4
from pdb import set_trace as bp

# from time import time as timer
# strt = timer()
class Reactions():
    def __init__(self,myGuy,Re,So,nReact,nComp,s_reactions,mySinks,mySources,areactants,RHS,order):
        self.myGuy = myGuy
        self.Re = Re
        self.So = So
        self.nReact = nReact
        self.nComp = nComp
        self.s_reactions = s_reactions
        self.mySinks = mySinks
        self.mySources = mySources
        self.areactants = areactants
        self.rhs = RHS
        self.nMaxR = 3
        self.nMaxP = 3
        self.groupsR = ['TOTAL']
        self.groupsP = ['TOTAL']
        self.order = order

    def lines(self):
        self.sources = []
        self.sinks = []
        self.sourceLabels = []
        self.sinkLabels = []
        groupsR = self.groupsR
        groupsP = self.groupsP
        nMaxR = min(self.nMaxR,sum(self.mySinks))
        nMaxP = min(self.nMaxP,sum(self.mySources))
        myGuy = self.myGuy
        # if self.myGuy in groupsR: groupsR.pop(groupsR.index(myGuy))
        myGuysBuddies = []
        for r in self.areactants[self.mySinks]:
            if all([x==myGuy for x in r]):
                myGuysBuddies.append(myGuy)
            else:
                for rr in r:
                    if rr != myGuy:
                        myGuysBuddies.append(rr)
        myGuysBuddies = np.array(myGuysBuddies)

        ###########################################################################################
        # Sinks / Reactivity
        ###########################################################################################
        if 'TOTAL' in groupsR:
            self.sinks.append(self.Re.sum(1))
            self.sinkLabels.append('Total React')
            groupsR.pop(groupsR.index('TOTAL'))
        if len(groupsR)>0 or nMaxR>0:
            for group in groupsR:
                indices = []
                if all([is_number(n) for n in group.split(',')]):
                    indices = [int(n) for n in group.split(',')]
                    indicesRed = reduceindices(indices, self.mySinks)
                    label = ' & '.join([self.s_reactions[i] for i in indices])
                else:
                    indicesRed = np.concatenate([np.where(myGuysBuddies == c)[0] for c in group.split(',')])
                    label = group
                if len(indicesRed)>0:
                    self.sinks.append(self.Re[:,indicesRed].sum(1).flatten())
                    self.sinkLabels.append(label)

            if len(groupsR)==0:
                sorter = np.argsort(self.Re.sum(0))[::-1]
                label = self.s_reactions[self.mySinks][sorter][:nMaxR]
                label = list(label) if nMaxR>1 else [label[0]]
                self.sinkLabels += label
                if nMaxR>1:
                    [self.sinks.append(y.flatten()) for y in self.Re[:,sorter][:,:nMaxR].T]
                else:
                    self.sinks.append(self.Re[:,sorter][:,:nMaxR].flatten())

        ###########################################################################################
        # Sources / Production
        ###########################################################################################
        if len(self.rhs[self.myGuy])==0:
            return
        if 'TOTAL' in groupsP:
            self.sources.append(self.So.sum(1))
            self.sourceLabels.append('Total production')
            groupsP.pop(groupsP.index('TOTAL'))

        if len(groupsP)>0 or nMaxP>0:
            for group in groupsP:
                indices=[]
                if all([is_number(n) for n in group.split(',')]):
                    indices = [int(n) for n in group.split(',')]
                    indicesRed = reduceindices(indices,self.mySources)
                    label = ' & '.join([self.s_reactions[i] for i in indices])
                    if len(indicesRed)>0:
                        self.sources.append(self.So[:,indicesRed].sum(1).flatten())
                        self.sourceLabels.append(label)

            if len(groupsP)==0:
                sorter = np.argsort(self.So.sum(0))[::-1]
                label = self.s_reactions[self.mySources][sorter][:nMaxP]
                label = list(label) if nMaxP>1 else [label[0]]
                self.sourceLabels += label
                if nMaxP>1:
                    [self.sources.append(y.flatten()) for y in self.So[:,sorter][:,:nMaxP].T]
                else:
                    self.sources.append(self.So[:,sorter][:,:nMaxP].flatten())
    # pass

def loadZip(f):
    file = gzip.GzipFile(f, 'rb')
    object = pickle.loads(file.read())
    file.close()
    return object

def reduceindices(iL,L,boolean=False):
    """Boolean vector and it's index mapped to Boolean vector"""
    a = np.array([True if i in iL else False for i in range(len(L))])[L]
    if boolean:
        return a
    else:
        return np.flatnonzero(a)

def getSto(Cin,time,myGuy,includeNull,chem,case,Dump=False):

    C = Cin.copy()

    def sourceR(name):
        i = dindices[name]
        if includeNull:
            return StoP[:,i]>0
        else:
            return Sto[:,i]>0
    def sinkR(name):
        i = dindices[name]
        if includeNull:
            return StoR[:,i]<0
        else:
            return Sto[:,i]<0
    def sourceSM(name):
        i = dindices[name]
        if includeNull:
            m = StoP[:,i]>0 # [sic!]
            return np.abs(np.where(StoR[m,:]<0,StoR[m,:],0))
        else:
            m = Sto[:,i]>0
            return np.abs(np.where(Sto[m,:]<0,Sto[m,:],0))
    def sinkSM(name):
        i = dindices[name]
        if includeNull:
            m = StoR[:,i]<0
            return np.abs(np.where(StoR[m,:]<0,StoR[m,:],0))
        else:
            m = Sto[:,i]<0
            return np.abs(np.where(Sto[m,:]<0,Sto[m,:],0))
    def sourceMu(name):
        i = dindices[name]
        if includeNull:
            m = StoP[:,i]>0
            return StoP[m,i]
        else:
            m = Sto[:,i]>0
            return Sto[m,i]
    def sinkMu(name):
        i = dindices[name]
        if includeNull:
            m = StoR[:,i]<0
            return np.abs(StoR[m,i])
        else:
            m = Sto[:,i]<0
            return np.abs(Sto[m,i])

    if Dump:
        process_reactions(chem=chem)

    nReact,nComp,Sto,dindices,s_reactants,s_reactions,areactants,StoR,StoP,RHS = loadZip(f'{chem}/Sto.zip')
    RC = np.genfromtxt(f'{case}/Reaction_rates.txt')
    nP, nR = sum(sourceR(myGuy)),sum(sinkR(myGuy))
    cMyGuy = np.ones(C.shape[0])*(C[:,dindices[myGuy]])

    homomolSecondOrder = np.arange(nReact)[[all([a==myGuy for a in ar]) for ar in areactants]]
    # bp()
    if len(homomolSecondOrder)>0:
        RC[:,homomolSecondOrder] = (RC[:,homomolSecondOrder].T*cMyGuy).T

    # Re = np.prod(np.transpose(np.tile(C[None,:,:] , (nR,1,1)), (1,2,0)),1, where=sinkSM(myGuy).T>0)*RC[:,sinkR(myGuy)]

    if includeNull:
        So = np.prod(np.transpose(np.tile(C[None,:,:] , (nP,1,1)), (1,2,0)),1, where=sourceSM(myGuy).T>0)*RC[:,sourceR(myGuy)]*sourceMu(myGuy)
        C[:,dindices[myGuy]] = 1.0
    else:
        C[:,dindices[myGuy]] = 1.0
        So = np.prod(np.transpose(np.tile(C[None,:,:] , (nP,1,1)), (1,2,0)),1, where=sourceSM(myGuy).T>0)*RC[:,sourceR(myGuy)]*sourceMu(myGuy)

    Re = np.prod(np.transpose(np.tile(C[None,:,:] , (nR,1,1)), (1,2,0)),1, where=sinkSM(myGuy).T>0)*RC[:,sinkR(myGuy)]*sinkMu(myGuy)

    return Reactions( myGuy, Re,So,nReact,nComp,s_reactions,sinkR(myGuy),sourceR(myGuy),areactants,RHS,np.sum(StoR,1) )

if __name__ == "__main__":
    from proc_reactivity import process_reactions,is_number
    chroot = '../../../src/chemistry/'
    ch = chroot+netCDF4.Dataset('/home/pecl/05-ARCA/ARCA-box/INOUT/STARTUP_0001/RKTV/Chemistry.nc').Chemistry_module
    myGuy = 'HO2'
    nMaxR,nMaxP = 5,5
    Dump = False
    includeNull=True
    C = loadZip('C.zip')
    nt = C.shape[0]
    time = np.linspace(0,12,nt)
    cs = '/home/pecl/05-ARCA/ARCA-box/INOUT/STARTUP_0001/RKTV'

    # Re,So,nReact,nComp,s_reactions,mySinks,mySources,areactants,RHS
    ROBJ = getSto(C,time, myGuy,includeNull,ch,cs,Dump)

    nMaxR = min(nMaxR,sum(ROBJ.mySinks))
    nMaxP = min(nMaxP,sum(ROBJ.mySources))
    # ProdReacFilter = [ True if i in ROBJ.RHS[myGuy] else False for i in range(ROBJ.nReact)]

    # groupsR = ['CO,CH4,SDD','O3','NO']
    groupsR = ['CO,CH4,O3','HO2','18']
    # groupsP = [1365, 'TOTAL']
    groupsP = ['35','TOTAL']

    groupsR = ['HO2']
    groupsP = ['TOTAL']

    import matplotlib.pyplot as plt
    plt.ion()

    f,(a1,a2) = plt.subplots(2)
    a1.set_title('Reactivity')
    a2.set_title('Production')
    # if myGuy in groupsR: groupsR.pop(groupsR.index(myGuy))

    myGuysBuddies = []
    for r in ROBJ.areactants[ROBJ.mySinks]:
        if all([x==myGuy for x in r]):
            myGuysBuddies.append(myGuy)
        else:
            for rr in r:
                if rr != myGuy:
                    myGuysBuddies.append(rr)

    myGuysBuddies = np.array(myGuysBuddies)

    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    # Sinks / Reactivity
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    if 'TOTAL' in groupsR:
        a1.plot(time, ROBJ.Re.sum(1), c='k', lw=1, ls='--', zorder=10,label='Total React')
        groupsR.pop(groupsR.index('TOTAL'))

    for group in groupsR:
        indices = []
        print(group.split(','))
        if all([is_number(n) for n in group.split(',')]):
            indices = [int(n) for n in group.split(',')]
            indicesRed = reduceindices(indices, ROBJ.mySinks)
            label = ' & '.join([ROBJ.s_reactions[i] for i in indices])
        else:
            indicesRed = np.concatenate([np.where(myGuysBuddies == c)[0] for c in group.split(',')])
            label = group
        if len(indicesRed)>0:
            a1.plot(time, ROBJ.Re[:,indicesRed].sum(1),label=label, lw=3)

    if len(groupsR)==0:
        sorter = np.argsort(ROBJ.Re.sum(0))[::-1]
        label = ROBJ.s_reactions[ROBJ.mySinks][sorter][:nMaxR]
        label = label if nMaxR>1 else label[0]
        a1.plot(time,ROBJ.Re[:,sorter][:,:nMaxR], label=label, lw=3)

    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    # Sources / Production
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    if len(ROBJ.rhs[myGuy])>0:
        if 'TOTAL' in groupsP:
            a2.plot(time, ROBJ.So.sum(1), c='k', lw=1, ls='--', zorder=10,label='Total Source')
            groupsP.pop(groupsP.index('TOTAL'))

        for group in groupsP:
            # if any([True if 'int' in str(type(g)) else False for g in groupsP]):
            if all([is_number(n) for n in group.split(',')]):
                indices = [int(n) for n in group.split(',')]
                indicesRed = reduceindices(indices,ROBJ.mySources)
                label = ' & '.join([ROBJ.s_reactions[i] for i in indices])
                if len(indicesRed)>0:
                    a2.plot(time,ROBJ.So[:,indicesRed].sum(1), label=label, lw=3)
        if len(groupsP)==0:
            sorter = np.argsort(ROBJ.So.sum(0))[::-1]
            label = ROBJ.s_reactions[ROBJ.mySources][sorter][:nMaxP]
            label = label if nMaxP>1 else label[0]
            a2.plot(time,ROBJ.So[:,sorter][:,:nMaxP], label=label, lw=3)

        a1.legend()
        a2.legend()
else:
    from modules.proc_reactivity import process_reactions,is_number
