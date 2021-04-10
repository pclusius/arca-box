#setup.py

import sys
import os

if sys.version_info.major < 3:
    print('Using ARCA Gui requires Python 3 (preferrably 3.6 or later)')
    quit()

if os.name == 'nt':
    le = '\r\n'
    posix = False
else:
    le = '\n'
    posix = True

import platform
operatingsystem = platform.system()
# "Windows"/"Linux"/"Darwin"
if '--csc' in sys.argv:
    csc = True
else:
    csc = False
out = 0
curr_path = os.path.split(os.getcwd())[0]
os.chdir(curr_path)

python = input('Type the command you use to call Python 3 (hit Enter for "python3"): ')
print(python)
if python == '':
    python = 'python3'


pyt = input('Install necessary Python packages? (y/n)?: ')

if pyt == 'y' or pyt == 'Y':
    print()
    print("Ok, let's see  what we need...")
    print()

    outpyt = os.system("%s -m pip install --user numpy scipy matplotlib requests"%python)
    outpyt = os.system("%s -m pip install --user netCDF4"%python)
    outpyt = os.system("%s -m pip install --user PyQt5"%python)
    outpyt = os.system("%s -m pip install --user pyqtgraph"%python)

    if outpyt != 0:
        upgr = input('Unfortunately the Python module installation did not work, updating setuptools could help.\nProceed and try again? (y/n)?: ')
        if upgr == 'y' or upgr == 'Y':
            outpyt = os.system("%s -m pip install --upgrade setuptools"%python)
            outpyt += os.system("%s -m pip install --user numpy scipy matplotlib netCDF4 PyQt5 pyqtgraph"%python)
        if outpyt != 0:
            print('Unfortunately still some Python modules failed to install.')

comp = input('Compile the Fortran module (y/n)?: ')

if comp == 'y' or comp == 'Y':
    print()
    print("This will take a while, have patience...")
    print()

    out = os.system("make")
    if out == 0:
        print("+------------------------------------------+\n")
        print("Compiling of arcabox.exe was succesful!\n")

if out == 0:
    if operatingsystem == 'Darwin':
        f = open('run_arca.command', 'w')
    else:
        f = open('run_arca.sh', 'w')

    f.write('#!/bin/bash%s'%le)
    if posix: f.write('cd %s %s'%(curr_path, le))
    f.write()
    if csc:
        f.write('module load python-data/3.7.6-1')
        f.write('module load gcc/9.1.0')
        f.write('module load netcdf/4.7.0')
        f.write('module load netcdf-fortran/4.4.4')
    else:
        print('ARCA user interface can be used on CSC with NoMachine. To get correct settings,')
        print('run this script with --csc flag or edit run_arca.sh and uncomment the csc options.')
        f.write('# If you are using the user interface on CSC, uncomment next 4 lines')
        f.write('# module load python-data/3.7.6-1')
        f.write('# module load gcc/9.1.0')
        f.write('# module load netcdf/4.7.0')
        f.write('# module load netcdf-fortran/4.4.4')
    f.write()
    f.write('%s ModelLib/gui/ARCA_gui.py %s'%(python, le))
    f.close()





    if posix:
            if operatingsystem == 'Darwin':
                out = os.system("chmod a+x run_arca.command")
                if out==0:
                    print("\nChanged the permission of 'run_arca.command'")
            else:
                out = os.system("chmod a+x run_arca.sh")
                if out==0:
                    print("\nChanged the permission of 'run_arca.sh'")

    if os.name == 'nt':
            print('\nSetup of GUI succesful. Run ARCA from Cygwin terminal with "sh run_arca.sh"')
    elif operatingsystem == "Linux":
        out = os.system("cp %s/ModelLib/gui/icons/thebox_ico.png ~/.icons/arca.png"%curr_path)
        if out==0:
            print("\nSuccesful copied ARCA icon to ~/.icons/")

        f = open('arca.desktop', 'w')
        f.write("""[Desktop Entry]
Name=ARCA Box Model
Comment=Atmospherically Relevant Chemistry and Aerosol Box model
Exec=%s/run_arca.sh
Icon=%s/ModelLib/gui/icons/thebox_ico.png
Terminal=true
Type=Application
Categories=GNOME;GTK;Core;
StartupWMClass=ARCAbox utility
"""%(curr_path,curr_path))
        f.close()
        print('\nSetup of GUI succesful. Run ARCA from terminal with "sh run_arca.sh", or by double clicking the ARCA box icon\n')
    elif operatingsystem == "Darwin":
        print('\nSetup of GUI succesful. Run ARCA from terminal with "./run_arca.command", or by double clicking the run_arca.command file.\n')

else:
    print("""Compiling the Fortran executable failed. First check that Fortran compiler is working.
The default compiler is gfortran, if some other compiler is used, change the variable "F90" in the
"makefile" to corresponding compiler. If the compiler is ok, make sure that netcdf-fortran and curl-dev
is installed from Cygwin. Email the developers if you keep having troubles.

You can rerun this script and omit the compiling of the Fortran model, this will complete the installation of thr GUI.
""")

#
