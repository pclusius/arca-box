import os,re,sys

changes = {'INDEX':'NUMBER'}
pattern = '|'.join(['^ *?%s'%k for k in changes.keys()])
n = len(changes)

def compatibility_layer(file):
    with open(file, 'r') as fi:
        with open(file+'.tmptmp', 'x') as fo:
            for line in fi:
                r = re.match(pattern, line)
                if r:
                    res = r.group().strip()
                    line = ' '+changes[res]+line[r.span()[1]:]

                    # for k in changes.keys():
                    #     line.replace(changes[k])
                fo.write(line)
    os.system("mv %s %s" %(file+'.tmptmp', file))
    # os.remove(file+'.tmptmp')
compatibility_layer(sys.argv[1])
