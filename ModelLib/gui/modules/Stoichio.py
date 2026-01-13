import numpy as np
import os, re, pickle, gzip,netCDF4
# from pdb import set_trace as bp
import psutil
# from memory_profiler import profile
# from time import time as timer
# strt = timer()
class Reactions():
    def __init__(self,myGuy,Re,So,nReact,nComp,s_reactions,mySinks,mySources,areactants,order):
        self.myGuy = myGuy
        self.Re = Re
        self.So = So
        self.nReact = nReact
        self.nComp = nComp
        self.s_reactions = s_reactions
        self.mySinks = mySinks
        self.mySources = mySources
        self.areactants = areactants
        # self.rhs = RHS
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
                    label = ' & '.join([f'R({i+1}): '+self.s_reactions[i] for i in indices])
                else:
                    indicesRed = np.concatenate([np.where(myGuysBuddies == c)[0] for c in group.split(',')])
                    label = group
                if len(indicesRed)>0:
                    self.sinks.append(self.Re[:,indicesRed].sum(1).flatten())
                    self.sinkLabels.append(label)

            if len(groupsR)==0:
                sorter = np.argsort(self.Re.sum(0))[::-1]
                label = [f'R({i+1}): '+s for i,s in zip(np.arange(self.nReact)[self.mySinks][sorter][:nMaxR],self.s_reactions[self.mySinks][sorter][:nMaxR])]
                label = list(label) if nMaxR>1 else [label[0]]
                self.sinkLabels += label
                if nMaxR>1:
                    [self.sinks.append(y.flatten()) for y in self.Re[:,sorter][:,:nMaxR].T]
                else:
                    self.sinks.append(self.Re[:,sorter][:,:nMaxR].flatten())

        ###########################################################################################
        # Sources / Production
        ###########################################################################################
        if self.mySources.sum()==0:
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
                    label = ' & '.join([f'R({i+1}): '+self.s_reactions[i] for i in indices])
                    if len(indicesRed)>0:
                        self.sources.append(self.So[:,indicesRed].sum(1).flatten())
                        self.sourceLabels.append(label)

            if len(groupsP)==0:
                sorter = np.argsort(self.So.sum(0))[::-1]
                label = [f'R({i+1}): '+s for i,s in zip(np.arange(self.nReact)[self.mySources][sorter][:nMaxP],self.s_reactions[self.mySources][sorter][:nMaxP])]
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

# @profile
def getSto(Cin,myGuy,includeNull,chem,case,Dump=False,include_homomolSecondOrder=True,jumpTime=1,loop=False):

    C = Cin.copy().astype(np.float32)

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
            return np.abs(np.where(StoR[m,:]<0,StoR[m,:],0)).astype(np.float32)
        else:
            m = Sto[:,i]>0
            return np.abs(np.where(Sto[m,:]<0,Sto[m,:],0)).astype(np.float32)
    def sinkSM(name):
        i = dindices[name]
        if includeNull:
            m = StoR[:,i]<0
            return np.abs(np.where(StoR[m,:]<0,StoR[m,:],0)).astype(np.float32)
        else:
            m = Sto[:,i]<0
            return np.abs(np.where(Sto[m,:]<0,Sto[m,:],0)).astype(np.float32)
    def sourceMu(name):
        i = dindices[name]
        if includeNull:
            m = StoP[:,i]>0
            return StoP[m,i].astype(np.float32)
        else:
            m = Sto[:,i]>0
            return Sto[m,i].astype(np.float32)
    def sinkMu(name):
        i = dindices[name]
        if includeNull:
            m = StoR[:,i]<0
            return np.abs(StoR[m,i]).astype(np.float32)
        else:
            m = Sto[:,i]<0
            return np.abs(Sto[m,i]).astype(np.float32)

    if Dump:
        process_reactions(chem=chem)

    nReact,nComp,Sto,dindices,s_reactions,areactants,StoR,StoP = loadZip(f'{chem}/Sto.zip')
    # del s_reactants

    RC = netCDF4.Dataset(f'{case}/Rates.nc').variables['Reaction_rates'][::jumpTime,:].data.astype(np.float32)

    sourceR_MG,sinkR_MG = sourceR(myGuy),sinkR(myGuy)
    nP, nR = sum(sourceR_MG),sum(sinkR_MG)
    GiB = C.shape[0]*C.shape[1]*max(nP,nR)*8 /2**30
    freeMem = psutil.virtual_memory().available /2**30
    if GiB / freeMem > 0.75:
    # if GiB>3:
        loop = True
        print(f'This requires too much memory ({GiB:0.2f} GiB > 75% of {freeMem:0.2f} GiB available), looping over time')
    if include_homomolSecondOrder:
        # cMyGuy = np.ones(C.shape[0])*(C[:,dindices[myGuy]])
        homomolSecondOrder = np.arange(nReact)[ [ar.count(myGuy)>1 for ar in areactants] ]
        if len(homomolSecondOrder)>0:
            RC[:,homomolSecondOrder] = (RC[:,homomolSecondOrder].T*C[:,dindices[myGuy]]).T
    # bp()
    RB = Reactions( myGuy, None,None,nReact,nComp,s_reactions,sinkR_MG,sourceR_MG,areactants,np.sum(StoR,1) )
    if not loop:
        if includeNull:
            RB.So = np.prod(np.transpose(np.tile(C[None,:,:] , (nP,1,1)), (1,2,0)),1, where=sourceSM(myGuy).T>0)*RC[:,sourceR_MG]*sourceMu(myGuy)
            C[:,dindices[myGuy]] = 1.0
        else:
            C[:,dindices[myGuy]] = 1.0
            RB.So = np.prod(np.transpose(np.tile(C[None,:,:] , (nP,1,1)), (1,2,0)),1, where=sourceSM(myGuy).T>0)*RC[:,sourceR_MG]*sourceMu(myGuy)

        RB.Re = np.prod(np.transpose(np.tile(C[None,:,:] , (nR,1,1)), (1,2,0)),1, where=sinkSM(myGuy).T>0)*RC[:,sinkR_MG]*sinkMu(myGuy)
    else:
        # print('Not enough memory, looping through time, this will take time')
        RB.So = np.zeros((C.shape[0],sourceR_MG.sum()))
        RB.Re = np.zeros((C.shape[0],sinkR_MG.sum()))
        print(f'Processing time indexes ({C.shape[0]} in total):')
        for it in range(C.shape[0]):
            if it%10==0:
                print(it, end='...',flush=True)
            if includeNull:
                RB.So[it,:] = np.sum(((C[it,:]* sourceSM(myGuy)).T * RC[it,sourceR_MG]) ,0)* sourceMu(myGuy)
                C[it,dindices[myGuy]] = 1.0
            else:
                C[it,dindices[myGuy]] = 1.0
                RB.So[it,:] = np.sum(((C[it,:]* sourceSM(myGuy)).T * RC[it,sourceR_MG]) ,0)* sourceMu(myGuy)
            RB.Re[it,:] = np.sum(((C[it,:]* sinkSM(myGuy)).T * RC[it,sinkR_MG]),0)*sinkMu(myGuy)
        print(' done.')
    # del C
    # del RC
    # del Sto,StoP
    return RB # Reactions( myGuy, Re,So,nReact,nComp,s_reactions,sinkR_MG,sourceR_MG,areactants,np.sum(StoR,1) )

if __name__ == "__main__":
    pass
    # from proc_reactivity import process_reactions,is_number, to_number
    #
    # chroot = '../../../src/chemistry/'
    # cs = '/home/pecl/05-ARCA/ARCA-box/INOUT/STARTUP_0001/INTTEST/'
    # cs = '/home/pecl/05-ARCA/ARCA-box/INOUT/HYDEAPRIL/PC_2018-04-11/NEWARCA2/'
    # nc = netCDF4.Dataset(f'{cs}/Chemistry.nc')
    # dimt = nc.dimensions['time'].size
    # if dimt>100:
    #     jumpTime = int(dimt/100)+1
    #
    # ch = chroot+nc.Chemistry_module
    # myGuy = 'OH'
    # nMaxR,nMaxP = 5,5
    # Dump = False
    # includeNull = True
    # saveMem = True
    # saveMem = False
    # # C = loadZip('C.zip')
    # C = nc.variables['CH_GAS'][::jumpTime,:].data
    # time = nc.variables['TIME_IN_HRS'][::jumpTime]
    # # time = np.linspace(0,12,nt)
    # nt = C.shape[0]
    #
    # # Re,So,nReact,nComp,s_reactions,mySinks,mySources,areactants,RHS
    # ROBJ = getSto(C,myGuy,includeNull,ch,cs,Dump=Dump,jumpTime=jumpTime,loop=saveMem)
    #
    # # ROBJ.nMaxR = min(nMaxR,sum(ROBJ.mySinks))
    # # ROBJ.nMaxP = min(nMaxP,sum(ROBJ.mySources))
    # # ProdReacFilter = [ True if i in ROBJ.RHS[myGuy] else False for i in range(ROBJ.nReact)]
    #
    # # groupsR = ['CO,CH4,SDD','O3','NO']
    # ROBJ.groupsR = ['CO,CH4,O3','HO2','18']
    # # groupsP = [1365, 'TOTAL']
    # ROBJ.groupsP = ['35','TOTAL']
    #
    # ROBJ.groupsR = ['TOTAL']
    # ROBJ.groupsP = ['TOTAL']
    #
    # import matplotlib.pyplot as plt
    # plt.ion()
    #
    # f,(a1,a2) = plt.subplots(2)
    # a1.set_title('Reactivity')
    # a2.set_title('Production')
    # ROBJ.lines()
    # for lab,lin in zip(ROBJ.sinkLabels,ROBJ.sinks):
    #     a1.plot(time, lin, label=lab)
    # for lab,lin in zip(ROBJ.sourceLabels,ROBJ.sources):
    #     a2.plot(time, lin, label=lab)
    # [a.legend() for a in [a1,a2]]
    #
    # ROBJ.groupsR = ['TOTAL']
    # ROBJ.groupsP = ['TOTAL']
    # ROBJ.lines()
    # iMyGuy = nc.variables['CH_GAS'].names.split(',').index(myGuy)
    # cc = C[:,iMyGuy]
    # import scipy as sc
    # f,(a1,a2) = plt.subplots(2)
    # for lab,lin in zip(ROBJ.sinkLabels,ROBJ.sinks):
    #     a1.plot(time[1:], sc.integrate.cumtrapz(lin*ROBJ.sources[0],time*3600), label=lab)
    # for lab,lin in zip(ROBJ.sourceLabels,ROBJ.sources):
    #     a2.plot(time[1:], sc.integrate.cumtrapz(lin,time*3600), label=lab)
    # [a.legend() for a in [a1,a2]]


    # sc.integrate.cumtrapz(y,x))

    # # if myGuy in groupsR: groupsR.pop(groupsR.index(myGuy))
    #
    # myGuysBuddies = []
    # for r in ROBJ.areactants[ROBJ.mySinks]:
    #     if all([x==myGuy for x in r]):
    #         myGuysBuddies.append(myGuy)
    #     else:
    #         for rr in r:
    #             if rr != myGuy:
    #                 myGuysBuddies.append(rr)
    #
    # myGuysBuddies = np.array(myGuysBuddies)
    #
    # ###########################################################################################
    # ###########################################################################################
    # ###########################################################################################
    # # Sinks / Reactivity
    # ###########################################################################################
    # ###########################################################################################
    # ###########################################################################################
    # if 'TOTAL' in groupsR:
    #     a1.plot(time, ROBJ.Re.sum(1), c='k', lw=1, ls='--', zorder=10,label='Total React')
    #     groupsR.pop(groupsR.index('TOTAL'))
    #
    # for group in groupsR:
    #     indices = []
    #     print(group.split(','))
    #     if all([is_number(n) for n in group.split(',')]):
    #         indices = [int(n) for n in group.split(',')]
    #         indicesRed = reduceindices(indices, ROBJ.mySinks)
    #         label = ' & '.join([f'R({i+1}): '+ROBJ.s_reactions[i] for i in indices])
    #     else:
    #         indicesRed = np.concatenate([np.where(myGuysBuddies == c)[0] for c in group.split(',')])
    #         label = group
    #     if len(indicesRed)>0:
    #         a1.plot(time, ROBJ.Re[:,indicesRed].sum(1),label=label, lw=3)
    #
    # if len(groupsR)==0:
    #     sorter = np.argsort(ROBJ.Re.sum(0))[::-1]
    #     label = [f'R({i+1}): '+s for i,s in zip(sorter[:nMaxR],ROBJ.s_reactions[ROBJ.mySinks][sorter][:nMaxR])]
    #     label = label if nMaxR>1 else label[0]
    #     a1.plot(time,ROBJ.Re[:,sorter][:,:nMaxR], label=label, lw=3)
    #
    # ###########################################################################################
    # ###########################################################################################
    # ###########################################################################################
    # # Sources / Production
    # ###########################################################################################
    # ###########################################################################################
    # ###########################################################################################
    # if len(ROBJ.rhs[myGuy])>0:
    #     if 'TOTAL' in groupsP:
    #         a2.plot(time, ROBJ.So.sum(1), c='k', lw=1, ls='--', zorder=10,label='Total Source')
    #         groupsP.pop(groupsP.index('TOTAL'))
    #
    #     for group in groupsP:
    #         # if any([True if 'int' in str(type(g)) else False for g in groupsP]):
    #         if all([is_number(n) for n in group.split(',')]):
    #             indices = [int(n) for n in group.split(',')]
    #             indicesRed = reduceindices(indices,ROBJ.mySources)
    #             label = ' & '.join([f'R({i+1}): '+ROBJ.s_reactions[i] for i in indices])
    #             if len(indicesRed)>0:
    #                 a2.plot(time,ROBJ.So[:,indicesRed].sum(1), label=label, lw=3)
    #     if len(groupsP)==0:
    #         sorter = np.argsort(ROBJ.So.sum(0))[::-1]
    #         label = [f'R({i+1}): '+s for i,s in zip(sorter[:nMaxP],ROBJ.s_reactions[ROBJ.mySources][sorter][:nMaxP])]
    #         label = label if nMaxP>1 else label[0]
    #         a2.plot(time,ROBJ.So[:,sorter][:,:nMaxP], label=label, lw=3)
    #
    #     a1.legend()
    #     a2.legend()
    # sc.integrate.cumtrapz(y,x))
else:
    from modules.proc_reactivity import process_reactions,is_number
