import re, sys,os

def updateINITFILE(file,updt=None):
    if updt is None:
        if os.system(f'grep "ARCA Box Model v." {file} 1>/dev/null') == 0:
            return 'File is already compatible with new MODS'
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        with open("../../required/version.txt", encoding="utf-8") as f:
            for line in f: break
        updt = 'ARCA Box Model v.'+line.replace('#','').strip('\n\r').strip('\n')

    # print(f'Checking that file {file} is compatible with current version')
    isOld = False
    os.system(f'mv {file} {file}_oldversion')
    counter = 1
    line_for_N_VARS = 0
    newfile_ = []
    with open(f'{file}_oldversion') as f:
        for line in f:
            a = re.search(r'MODS\(\d+\)', line)
            b = re.search(r'# +Created at:', line)
            c = re.search(r'NAMESDAT =', line)
            if a is not None and isOld:
                line_ = list(re.findall(r'(MODS\()(\d+)(\))(.+)(!)(.+)', line)[0])
                if '!' not in line_:
                    isOld = False
                    newfile_.append(line)
                    continue
                iex = line_.index('!')
                line_[iex+1] = "'"+line_[iex+1].strip()+"'"
                line_.pop(iex)
                line_[1] = f'{counter}'
                newfile_.append(''.join(line_)+'\n')
                counter += 1
            elif b is not None:
                newfile_.append(f'# {updt}\n')
                newfile_.append(line.replace('Created at: ','Created on '))
                isOld = True
            else:
                newfile_.append(line)
            if c is not None:
                line_for_N_VARS = len(newfile_)
                newfile_.append(" N_VARS = ")
    newfile_[line_for_N_VARS] += f'{counter-1}\n'
    newfile = open(file, 'w+')
    newfile.write(''.join(newfile_))
    newfile.close()

    return f'Edited {file} to work with current version, old file is renamed "{file}_oldversion"'


def parse_chamistry_NAMES(path):
    dir,file = os.path.split(path)
    names = ['TEMPK','PRESSURE','REL_HUMIDITY','CONDENS_SINK','CON_SIN_NITR',
        'SW_RADIATION','ION_PROD_RATE','NUC_RATE_IN','# DONT_CHANGE Anything ABOVE THIS LINE']
    with open(os.path.join(dir,'second_Parameters.f90'), 'r') as params:
        for line in params:
            m = re.findall('(INTEGER, PARAMETER :: ind)(_)([a-zA-Z0-9_]+)', line.strip())
            if len(m)>0:
                name = m[0][-1]
                if name.upper() != 'DUMMY':
                    names.append(name)
    if len(names)>0:
        with open(path, 'w+') as NAMES:
            [NAMES.write(name+'\n') for name in names]
    return f'Parsed new NAMES.dat of the compounds from {dir}'

def convertMD(file):
    empty = ''
    with open(file,'r') as f:
        for line in f:
            ln = line.strip()+'¤'
            empty += re.sub('^[-*]','¤*',ln)
    text = re.sub('(#+? )', r'££££££\1',empty)
    text = re.sub('¤¤+','££££££',text)
    text = re.sub('¤',' ',text)
    text = re.sub('££££££','¤¤',text)
    text = text.split('¤¤')
    retarr = []
    inList = False
    for i_t,t in enumerate(text):
        if inList and re.search(r'^\* ',t) is None:
            retarr.append('</ul>')
            inList = False
        if re.search(r'^\* ',t):
            t = re.sub(r'^\*', '<li>', t) + '</li>'
            if not inList:
                retarr.append('<ul>')
                inList = True
        elif re.search('^# ', t):
            t = re.sub('^# ', '<h1>', t)
            t += '</h1>'
        elif re.search('^## ', t):
            t = re.sub('^## ', '<h2>', t)
            t += '</h2>'
        elif re.search('^### ', t):
            t = re.sub('^### ', '<h3>', t)
            t += '</h3>'
        elif re.search('^#### ', t):
            t = re.sub('^#### ', '<h4>', t)
            t += '</h4>'
        else:
            t = '<p>' + t + '</p>'
        retarr.append(t)
    return ' '.join(retarr)

if __name__ == '__main__':
    print("""
    functions:
        convertMD(file)
        parse_chamistry_NAMES(path)
        updateINITFILE(file,updt=None)
    """)
