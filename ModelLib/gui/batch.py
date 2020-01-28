from os import path, listdir
from datetime import datetime, timedelta

def paths(addr):
    if path.exists(addr):
        if path.isdir(addr):
            return 0
        else:
            return 1
    else:
        return -1

def batch(begin = '2020-01-01', end='2020-01-05', case='Melpitz', run='X1', common_root=''):
    if common_root != '':
        if common_root[-1] != '/':
            common_root = common_root+'/'
        if paths(common_root[:-1]) == -1:
            return 'Common root does not exist. For safety this directory must exist beforehand.'

    if 'int' in str(type(begin)):
        if begin > end:
            return 'Ending index cannot be smaller than beginning'
        ind_range = range(begin,end+1)
        fmt = "%0"+str(max(4,len(str(end))))+"d"
        index_strings =[fmt %i for i in ind_range]
    if 'str' in str(type(begin)):
        begin = datetime.strptime(begin, "%Y-%m-%d")
        end = datetime.strptime(end, "%Y-%m-%d")
        if begin > end:
            return 'End date cannot be before beginning'

        ind_range = range((end-begin).days+1)
        index_strings = [datetime.strftime(begin+timedelta(days=i), "%Y-%m-%d") for i in ind_range]

    case = case.upper()
    run = run.upper()


    dirs_to_create = []
    files_to_create = []
    files_to_overwrite = []
    existing_runs = []
    conflicting_names = []

    for i in index_strings:

        # indir = '%s%s_%s/input_%s_%s' %(common_root,case,i,case,i)
        indir = '%s%s_%s/input' %(common_root,case,i)
        outdir = '%s%s_%s' %(common_root,case,i)
        inifile = indir+'/%s_%s_%s.conf' %(case,i,run)
        files_to_create.append(inifile)
        if paths(inifile) == 1:
            files_to_overwrite.append(inifile)

        for j,p in enumerate((outdir,indir, outdir+'/'+run)):
            if paths(p) < 0: dirs_to_create.append(p)
            elif paths(p) > 0 : conflicting_names.append(p)
            else:
                if j==2 and listdir(p) != []:
                    existing_runs.append(p)

    return dirs_to_create, conflicting_names, files_to_create, files_to_overwrite, existing_runs, index_strings, outdir+'/'+run
# batch()
def tagparser(tag, index):
    tag = tag.replace('<','')
    tag = tag.replace('>','')
    n = len(tag)
    if n>0:
        c = tag[0]
        try:
            index = int(index)
        except: pass
        if all(tag[i]==c for i in range(n)):
            if 'str' in str(type(index)):
                if c.upper() == 'Y' and n<5:
                    retv = index[4-n:4]
                elif c.upper() == 'M' and n<3:
                    retv = index[7-n:7]
                elif c.upper() == 'D' and n<3:
                    retv = index[10-n:]
                else:return '<unsupported tag for date>'
            elif 'int' in str(type(index)) and c.upper() == 'I':
                retv = ('%%0%dd'%n) %index
            else:return '<unsupported tag for number>'
            return retv
        else:return '<unsupported form>'
    else: return '<>'
