import numpy as np
import pickle, gzip,netCDF4

def is_number(s):
    try: return int(s)
    except: return None
def saveZip(o,f):
    file = gzip.GzipFile(f, 'wb')
    file.write(pickle.dumps(o, -1))
    file.close()
def loadZip(f):
    file = gzip.GzipFile(f, 'rb')
    object = pickle.loads(file.read())
    file.close()
    return object

def process_reactions(chem,ncfile=None):
    dnames,anames = {},[]
    dindices,aindices = {},[]
    dreactions,s_reactions = {},[]
    dreactants,s_reactants = {},[]
    dproducts,s_products = {},[]
    nComp = 0
    with open(f'{chem}/SPC_NAMES.txt', 'r') as file:
        for line in file:
            dindices[line.strip()] = nComp
            dnames[nComp] = line.strip()
            nComp += 1

    nReact = 0
    with open(f'{chem}/EQN_NAMES.txt', 'r') as file:
        for line in file:
            reaction = line.strip()
            left,right = reaction.split('-->')
            left = left.replace('+', '').strip().split()
            right = right.replace('+', '').strip().split()
            for il,c in enumerate(left):
                N = is_number(left[il])
                if N:
                    left[il] = left[il+1]
                    _ = [left.append(left[il+1]) for x in range(N-2)]
            for il,c in enumerate(right):
                N = is_number(right[il])
                if N:
                    right[il] = right[il+1]
                    _ = [right.append(right[il+1]) for x in range(N-2)]
            dreactions[nReact] = reaction
            dreactants[nReact] = left
            dproducts[nReact] = right
            s_reactions.append(reaction)
            s_reactants.append(' + '.join(left))
            s_products.append(' + '.join(right))
            nReact += 1

    s_reactions = np.array(s_reactions)
    s_reactants = np.array(s_reactants)
    s_products = np.array(s_products)

    Sto = np.zeros((nReact,nComp))
    StoR = np.zeros((nReact,nComp))
    StoP = np.zeros((nReact,nComp))

    for r in range(nReact):
        for c in dreactants[r]:
            Sto[r, dindices[c]] -= 1
            StoR[r, dindices[c]] -= 1
        for c in dproducts[r]:
            Sto[r, dindices[c]] += 1
            StoP[r, dindices[c]] += 1

    if ncfile is not None:
        nc = netCDF4.Dataset( ncfile )
        t = nc.variables['TIME_IN_HRS'][:]
        nt = len(t)
        C = np.zeros((nt,nComp))
        for i in range(nComp):
            C[:,i] = nc.variables[dnames[i]][:]

    areactants = np.array([v for v in dreactants.values()], dtype=object)
    aproducts = np.array([v for v in dproducts.values()], dtype=object)

    RHS = {n:[] for n in dindices.keys() }
    for k in dproducts.keys():
        for n in list(set(dproducts[k])):
            RHS[n].append(k)

    # LHS = {n:[] for n in dindices.keys() }
    # for k in dproducts.keys():
    #     for n in list(set(dreactants[k])):
    #         LHS[n].append(k)

    saveZip([nReact,nComp,Sto,dindices,s_reactants,s_reactions,areactants,StoR,StoP,RHS],f'{chem}/Sto.zip')

    if ncfile is not None:
        saveZip(C,'C.zip')
