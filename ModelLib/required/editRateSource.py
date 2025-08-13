import sys,os,re,shutil

f = sys.argv[1]
unchanged = True
if int(sys.version_info.major) < 3:
    exit('this script only works with Python3')

hits = 0
ratesfile = ''
with open(os.path.join(f,'second_Rates.f90'), 'r') as file:
    for line in file:
        if re.search(r"RCONST\s*=\s*RCONST\s*\*\s*R_F",line.strip('\n').upper()):
            hits +=1
        if re.search(r'^\s*END\s*SUBROUTINE\s*UPDATE_RCONST', line.upper()) and hits == 0:
            ratesfile += '  RCONST = RCONST * R_F\n'
        ratesfile += line
    if hits>1:
        shutil.move(os.path.join(f,'second_Rates.f90'),os.path.join(f,'second_Rates_CHECK.f90'))
        print(''.join(['\n']*10)+'    second_Rates.f90 is messed up, check manually and remove unnecessary RCONST manipulations')
        exit('    second_Rates.f90 is renamed to second_Rates_CHECK.f90 in %s'%f+''.join(['\n']*10))

if hits==0:
    mod_Rfile = open(os.path.join(f,'second_Rates.f90'), "w")
    mod_Rfile.write(ratesfile)
    mod_Rfile.close()
    unchanged = False
    print(' Edited second_Rates.f90 to accommodate manipulation of reaction constants from outside')

hits = 0
mainfile = ''
with open(os.path.join(f,'second_Main.f90'), 'r') as file:
    for line in file:
        if re.search(r"RCONST\s*=\s*RCONST\s*\*\s*R_F",line.strip('\n').upper()):
            hits +=1
        elif re.search(r"REAL\(DP\),\s*SAVE\s*::\s*R_F\(NREACT\)",line.strip('\n').upper()):
            hits +=1
        elif re.search(r"\s*!\s*MODIFY\s*RATE\s*CONSTANTS",line.strip('\n').upper()):
            hits +=1
        elif re.search(r'SUBROUTINE KPP_SETUP\(\)', line.strip('\n').upper()) != None:
            mainfile += '  SUBROUTINE KPP_SetUp(R_F_in)\n  REAL(DP), OPTIONAL :: R_F_in(NREACT)\n  IF (PRESENT(R_F_in)) R_F = R_F_in\n'
            hits +=1

        else:
            mainfile += line

if hits>0:
    mod_Rfile = open(os.path.join(f,'second_Main.f90'), "w")
    mod_Rfile.write(mainfile)
    mod_Rfile.close()
    unchanged = False
    print('Fixed old procedures in second_Main.f90 to v 1.3.2 (& >>)')


hits = 0
mainfile = ''
search_on = True
with open(os.path.join(f,'second_Global.f90'), 'r') as file:
    for line in file:
        if re.search(r"REAL\(DP\),\s*SAVE\s*::\s*R_F\(NREACT\)\s*=\s*1D0",line.strip('\n').upper()):
            hits +=1

        if re.search(r'!\s*INLINED\s*GLOBAL\s*VARIABLE\s*DECLARATIONS', line.upper()) and hits == 0 and search_on:
            mainfile += '  ! Enables changing reactions rates from the outside of chemistry\n'
            mainfile += '  REAL(DP), SAVE :: R_F(NREACT) = 1D0\n'
            search_on = False

        mainfile += line

if hits==0:
    mod_Rfile = open(os.path.join(f,'second_Global.f90'), "w")
    mod_Rfile.write(mainfile)
    mod_Rfile.close()
    unchanged = False
    print(' Edited second_Global.f90 to accommodate manipulation of reaction constants from outside')

if hits>1:
    shutil.move(os.path.join(f,'second_Global.f90'),os.path.join(f,'second_Global_CHECK.f90'))
    print(''.join(['\n']*10)+'    second_Global.f90 is messed up, check manually and remove unnecessary R_F declarations')
    exit('    second_Global.f90 is renamed to second_Global_CHECK.f90 in %s'%f+''.join(['\n']*10))

if unchanged:
    print('  ---')
