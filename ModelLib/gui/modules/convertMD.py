import re

def convertMD(text):
    for i_t,t in enumerate(text):
        # print(t)
        inList = False
        if re.search(r'^\* ',t):
            print('*')
            t = re.sub(r'^\*', '<li>', t)
            if not inList:
                t = '<ul>'+t
                inList = True
        elif inList:
            print('*inlis')
            text[i_t-1] += text[i_t-1]+'</ul>'
            inList = False
        if re.search('^# ', t):
            print('h1')
            t = re.sub('^# ', '<h1>', t)
            t += '</h1>'
        elif re.search('^## ', t):
            print('h2')
            text[i_t-1] += '</p>'
            t = re.sub('^## ', '<h2>', t)
            t += '</h2>'
        elif re.search('^### ', t):
            print('h3')
            t = re.sub('^### ', '<h3>', t)
            t += '</h3>'
        elif re.search('^#### ', t):
            print('h4')
            t = re.sub('^#### ', '<h4>', t)
            t += '</h4>'
        else:
            print('p')
            t = '<p>' + t + '</p4>'
        text[i_t] = t
    return '\ņ'.join(text)

empty = ''
with open("../../../ModelLib/required/version.txt",'r') as f:
    for line in f:
        ln = line.strip()+'¤'
        empty += re.sub('^[-*]','¤*',ln)
t = re.sub('¤¤+','££££££',empty)
t = re.sub('¤',' ',t)
t = re.sub('££££££','¤¤',t)
t = t.split('¤¤')
# x = convertMD(Y)
# print()
# f.close()

# t = '\n'.join(convertMD(text))

# t = re.sub('\ņ ?\n','¤', tx)

# print(t.replace('¤','\ņ'))
