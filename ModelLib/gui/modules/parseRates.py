import os, sys, re

path = ''
fact_str = '1d0'
def es(x, p=1):
    st = '%0.1e'%(x)
    if p==1:
        print(st)
    else:
        return st

def parseF():
    try:
        float(factor)
        fact_str = es(float(factor), 0).replace('e','d')
    except:
        fact_str = factor.strip()
    return fact_str

if len(sys.argv)>1:
    path = sys.argv[1]
    if len(sys.argv)>2:
        factor = sys.argv[2]
        fact_str = parseF()
else:
    print('Give the path to the chemistry directory as cmdline option. the script creates a namelist of the file.')
    path = input('You can also type the path here: \n')
    if len(path.split())>1:
        path, factor = path.split()
        fact_str = parseF()

if path == '':
    exit('No path given, exiting.')

elif not os.path.exists(path):
    exit('The given path does not exist. Exiting')

elif not os.path.exists(os.path.join(path, 'second_Rates.f90')):
    exit('The given path exists but does not contain a chemistry scheme. Exiting')


rates = open(os.path.join(path, 'RATES.dat'), 'w+')
ratesList = []
# rates.write('&NML_RATES\n')
# ratesList.append('&NML_RATES\n')
inrates = False
with open(os.path.join(path, 'second_Rates.f90')) as file:
    for ln in file:
        if re.search('SUBROUTINE Update_RCONST',ln) != None and inrates == False:
            inrates=True
        if re.search('END SUBROUTINE',ln) != None and inrates == True:
            break
        if re.search(r'RCONST\(\d+\) = \(.*\)', ln.strip('\n')) != None:
            a, b = ln.strip('\n').split(' = ')
            # rates.write(a + ' = ' + fact_str + ' ! * '+ b + '\n')
            ratesList.append(a + ' = ' + fact_str + ' ! * '+ b + '\n')
        elif re.search(r'\!( *)RCONST\(\d+\) = ', ln.strip('\n')) != None:
            # rates.write(ln)
            ratesList.append(ln)
# rates.write('/\n')
# ratesList.append('/\n')
ratesList = list(set(ratesList))
sortedRates = ['\n']*len(ratesList)
for s in ratesList:
    m = re.search(r'RCONST\((\d+)\)',s)
    if m:
        sortedRates[int(m.group(1))-1] = s

ratesList = []
inrates = False
with open(os.path.join(path, 'second_Initialize.f90')) as file:
    for ln in file:
        if re.search('SUBROUTINE Initialize',ln) != None and inrates == False:
            inrates=True
        if re.search('END SUBROUTINE',ln) != None and inrates == True:
            break
        if re.search(r'RCONST\(\d+\) = .*', ln.strip('\n')) != None:
            a, b = ln.strip('\n').split(' = ')
            # rates.write(a + ' = ' + fact_str + ' ! * '+ b + '\n')
            ratesList.append('! '+a.strip() + ' = ' + fact_str + ' ! * '+ b + '\n')
        elif re.search(r'\!( *)RCONST\(\d+\) = ', ln.strip('\n')) != None:
            # rates.write(ln)
            ratesList.append(ln)

ratesList = list(set(ratesList))
for s in ratesList:
    m = re.search(r'RCONST\((\d+)\)',s)
    if m:
        sortedRates[int(m.group(1))-1] = s

rates.write('&NML_RATES\n'+''.join(sortedRates)+'/\n')
rates.close()

print('Saved RATES.dat in ' + path)


#
