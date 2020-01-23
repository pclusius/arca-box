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
            return 'Common root does not exist. For safety\nthis directory must exist beforehand.'

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

    # roots = [case]#, case+'/Input']

    # for r in roots:
    #     r = common_root + r
    #     if paths(r)<0:
    #         dirs_to_create.append(r)
    #     elif paths(r)>0:
    #         conflicting_names.append(r)


    # FORMAT EXAMPLES
    # indir = 'common_root/case_xxxx-xx-xx/Input_case_xxxx-xx-xx'
    # outdir = 'common_root/case_xxxx-xx-xx'
    # rundir = 'common_root/case_xxxx-xx-xx/run'
    # inifile = 'common_root/case_xxxx-xx-xx/run/case_xxxx-xx-xx_run.conf'
    for i in index_strings:
        # inifile = inifile_f %(i,i)
        # indir = indir_f %(i)
        # outdir = outdir_f %(i)

        indir = '%s%s_%s/Input_%s_%s' %(common_root,case,i,case,i)
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

    return dirs_to_create, conflicting_names, files_to_create, files_to_overwrite, existing_runs, index_strings
# batch()
