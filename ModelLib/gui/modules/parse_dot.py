import os, sys
dot = []
filter=''
dir = 'home'
plot = True

if len(sys.argv)>1:
    dir = sys.argv[1]
if len(sys.argv)>2:
    filter= sys.argv[2]

os.system(f'cat {dir}/second_Monitor.f90 |grep -F -- "-->" >{dir}/chem.txt')
with open(f'{dir}/chem.txt','r') as file:
    for il,line in enumerate(file):
        if filter!=''and filter not in line:
            continue
        Re = f'R{il+1}'
        # dot.append(f'{Re} [class="reaction" shape=box style=filled fillcolor=lightblue]')
        left,right = line.strip().split('-->')
        left = left.replace("'","").replace(' 2 ',' ').replace(' 3 ',' ').replace(' 4 ',' ')
        right = right.split("'")[0].replace(' 2 ',' ').replace(' 3 ',' ').replace(' 4 ',' ')
        if right.strip()=='':right='dummy'

        # stl = '{'+' '.join([s.strip()+';' for s in left.split('+')])+'}'+f'-> {Re}'
        # dot.append(f'{stl.replace(";}","}")}')
        # str = f'{Re} -> '+'{'+' '.join([s.strip()+';' for s in right.split('+')])+'}'
        # dot.append(f'{str.replace(";}","}")}')

        # dot.append(f'{st} [color=blue]')
        # dot.append(f'{st} [color=blue]')

        # for l in left.split('+'):
        #     dot.append(f'{l.strip()} -> {Re} [color=red]')
        # for r in right.split('+'):
        #     dot.append(f'{Re} -> {r.strip()} [color=blue]')
        for l in left.split('+'):
            # dot.append(f'{l.strip()} -> {Re} [color=red]')
            for r in right.split('+'):
                dot.append(f'{l.strip()} -> {r.strip()}')


with open(f'{dir}/x.dot','w') as file:
    file.write('digraph chem {\n')
    # file.write('    layout="sfdp";\n')
    # # file.write('    beautify=true;\n')
    # file.write('    overlap="scale";\n')
    # file.write('    splines = true;\n')
    for line in dot:
        file.write(line+' ; \n')
    file.write('}\n')

if plot:
    os.system(f"dot -Tpng -Granksep=0.5 {dir}/x.dot >{dir}/x.png")
    # os.system("dot  -Tpng -Granksep=2 -Gnodesep=1 -Grankdir=LR -Gsplines=ortho -Nshape=box x.dot >x.png")
    os.system(f'open {dir}/x.png')
