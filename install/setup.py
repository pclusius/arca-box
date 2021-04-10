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
    outpyt = os.system("module load gcc/9.1.0")
    outpyt = outpyt + os.system("module load netcdf/4.7.0")
    outpyt = outpyt + os.system("module load netcdf-fortran/4.4.4")
    if outpyt != 0: print("Unable to load necessary modules, if you are not on CSC, this might not be a problem.")
else:
    csc = False
out = 0
curr_path = os.path.split(os.getcwd())[0]
os.chdir(curr_path)

python = input('Type the command you use to call Python 3 (hit Enter for "python3"): ')
print(python)
if python == '':
    python = 'python3'

if not os.path.exists('makefile'):
    # write makefile
    f = open('makefile', 'w')
    f.write("""
# compiler
F90 = gfortran

COPTI = -O1
BOPTI = -O3
PROF = -pg
PROF =
# Put .o and .mod files here:
 OBJDIR = build
 SRCDIR = src
 CHMDIR = Hyde2

$(shell mkdir -p $(OBJDIR)/$(CHMDIR))
# When compiling, search for files in these directories:
VPATH = $(OBJDIR):src:src/ACDC/ACDC_module_2016_09_23:src/ACDC/ACDC_module_ions_2018_08_31:Aerosol:$(OBJDIR)/$(CHMDIR)
FILE=ModelLib/required/version.txt
VERSION=`head -1 $(FILE)`

# Options reminders:
# -w suppresses warning messages

BOX_OPTS = -g -Wno-unused $(BOPTI) -ffree-line-length-none -cpp -DLINUX -DCHEM=\\"$(CHMDIR)\\" -DVERSION=\\"$(VERSION)\\" -DISACDC -J$(OBJDIR) -I$(OBJDIR)/$(CHMDIR) -I$(OBJDIR) -fcheck=bounds,do -Wall -Wextra -Wsurprising \\
-Wno-unused-dummy-argument -Wno-maybe-uninitialized -Wtabs -Wno-tabs -Wno-character-truncation -fbacktrace -ffpe-trap=invalid,zero,overflow $(PROF) -g -fcheck=all

CHEM_OPTS = -w -cpp $(PROF) $(COPTI) -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -J$(OBJDIR)/$(CHMDIR) -I$(OBJDIR)/$(CHMDIR)

ACDC_OPTS = -ffree-line-length-none -cpp -J$(OBJDIR) -I$(OBJDIR) -fcheck=all -ffpe-trap=invalid,zero,overflow -O3

CHEM_OBJECTS = $(addprefix $(OBJDIR)/$(CHMDIR)/, second_Precision.o second_Parameters.o second_Initialize.o second_Util.o second_Monitor.o second_JacobianSP.o \\
               second_LinearAlgebra.o second_Jacobian.o second_Global.o second_Rates.o second_Integrator.o second_Function.o \\
               second_Model.o second_Main.o second_reactivity.o)

BOX_OBJECTS = $(addprefix $(OBJDIR)/, constants.o auxillaries.o input.o solve_bases.o chemistry.o psd_scheme.o aerosol_dynamics.o output.o custom_functions.o)

PSD_OBJECTS = $(addprefix $(OBJDIR)/, constants.o input.o chemistry.o)

AEROSOL_OBJECTS = $(addprefix $(OBJDIR)/, constants.o input.o auxillaries.o)

ACDC_OBJECTS = $(addprefix $(OBJDIR)/, vodea.o vode.o acdc_system_AN_ions.o acdc_system_extras.o monomer_settings_acdc_NH3_ions.o solution_settings.o driver_acdc_J_ions.o \\
             acdc_equations_AN_ions.o get_acdc_J_ions.o)

ACDC_D_OBJECTS = $(addprefix $(OBJDIR)/, vodea.o vode.o acdc_system_AD_new.o monomer_settings_acdc_DMA.o solution_settings.o driver_acdc_D.o \\
                  acdc_equations_AD_new.o get_acdc_D.o)

# ---------------------------------- NETCDF configuring ---------------------------------------
# To get Netcdf working you need to use system dependent settings
# If you get errors related to LIBNET, this is the place to start troubleshooting

# This seems to be default on Linux Mint/Ubuntu laptop and Windows with Cygwin%s"""%le)
    if (operatingsystem == 'Linux' or operatingsystem == 'Windows') and not csc:
        f.write('NETLIBS =  -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdf  -lnetcdff -lcurl%s'%le)
    else:
        f.write('# NETLIBS =  -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdf  -lnetcdff -lcurl%s'%le)

    f.write('%s'%le)
    f.write('# On Mac, this might work ---------------------------------------%s'%le)
    if operatingsystem == 'Darvin':
        f.write('NETLIBS = -I/opt/local/include -L/opt/local/lib -L/usr/lib -lnetcdff -lnetcdf -lcurl -lhdf5 -lhdf5_hl%s'%le)
    else:
        f.write('# NETLIBS = -I/opt/local/include -L/opt/local/lib -L/usr/lib -lnetcdff -lnetcdf -lcurl -lhdf5 -lhdf5_hl%s'%le)

    f.write('%s'%le)
    f.write('# Or this%s'%le)
    f.write('#NETLIBS = -I$(NETCDF_INCLUDE) -L$(NETCDF_LIB) -L$(H5_LIB) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl%s'%le)


    f.write('%s'%le)
    f.write('# On CSC/Puhti ------------------------------------------------------%s'%le)
    if csc:
        f.write('NETLIBS =  -I/appl/spack/install-tree/gcc-9.1.0/netcdf-fortran-4.4.4-4tj6tj/include -L/appl/spack/install-tree/gcc-9.1.0/netcdf-fortran-4.4.4-4tj6tj/lib -lnetcdff -lnetcdf -lnetcdf%s'%le)
    else:
        f.write('# NETLIBS =  -I/appl/spack/install-tree/gcc-9.1.0/netcdf-fortran-4.4.4-4tj6tj/include -L/appl/spack/install-tree/gcc-9.1.0/netcdf-fortran-4.4.4-4tj6tj/lib -lnetcdff -lnetcdf -lnetcdf%s'%le)

    f.write("""
# Load the correct modules before compiling and running:
# module load python-data/3.7.6-1
# module load gcc/9.1.0
# module load netcdf/4.7.0
# module load netcdf-fortran/4.4.4
# ------------------------------- end NETCDF configuring ---------------------------------------

all: arcabox.exe

# Here is the link step:
arcabox.exe: arcabox.o $(BOX_OBJECTS) $(CHEM_OBJECTS) $(ACDC_OBJECTS) $(ACDC_D_OBJECTS) $(PSD_OBJECTS)
	$(F90) $(BOX_OPTS) $^ -o $@ $(NETLIBS)


# Here are the compile steps
# Main program
$(OBJDIR)/arcabox.o: src/ARCA_main.f90 $(CHEM_OBJECTS) $(BOX_OBJECTS) $(ACDC_OBJECTS) $(ACDC_D_OBJECTS) $(PSD_OBJECTS)
	 $(F90) $(BOX_OPTS) -c $< -o $@

#PSD_scheme
$(OBJDIR)/psd_scheme.o: src/PSD_scheme.f90 $(PSD_OBJECTS)
	$(F90) $(BOX_OPTS) -c $< -o $@

#Aerosol dynamic
$(OBJDIR)/aerosol_dynamics.o: $(SRCDIR)/aerosol_dynamics.f90 $(AEROSOL_OBJECTS)
	$(F90) $(BOX_OPTS) -c $< -o $@

$(OBJDIR)/solve_bases.o: src/solve_bases.f90 $(ACDC_OBJECTS) $(ACDC_D_OBJECTS)
	 $(F90) $(BOX_OPTS) -c $< -o $@

# Chemistry
$(OBJDIR)/$(CHMDIR)/second_Precision.o: chemistry/$(CHMDIR)/second_Precision.f90
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Parameters.o: chemistry/$(CHMDIR)/second_Parameters.f90 second_Precision.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Global.o: chemistry/$(CHMDIR)/second_Global.f90 second_Parameters.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Initialize.o: chemistry/$(CHMDIR)/second_Initialize.f90 second_Parameters.o second_Global.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_reactivity.o: chemistry/$(CHMDIR)/second_reactivity.f90 second_Precision.o second_Parameters.o second_Global.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Util.o: chemistry/$(CHMDIR)/second_Util.f90 second_Parameters.o second_Monitor.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Monitor.o: chemistry/$(CHMDIR)/second_Monitor.f90
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_JacobianSP.o: chemistry/$(CHMDIR)/second_JacobianSP.f90
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_LinearAlgebra.o: chemistry/$(CHMDIR)/second_LinearAlgebra.f90 second_Parameters.o second_JacobianSP.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Jacobian.o: chemistry/$(CHMDIR)/second_Jacobian.f90 second_Parameters.o second_JacobianSP.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Rates.o: chemistry/$(CHMDIR)/second_Rates.f90 second_Parameters.o second_Global.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Integrator.o: chemistry/$(CHMDIR)/second_Integrator.f90 second_Precision.o second_Global.o second_Parameters.o second_JacobianSP.o \\
second_LinearAlgebra.o second_Function.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Function.o: chemistry/$(CHMDIR)/second_Function.f90 second_Parameters.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Model.o: chemistry/$(CHMDIR)/second_Model.f90 second_Precision.o second_Parameters.o second_Global.o second_Function.o \\
second_Integrator.o second_Rates.o second_Jacobian.o second_LinearAlgebra.o second_Monitor.o second_Util.o second_Rates.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@
$(OBJDIR)/$(CHMDIR)/second_Main.o: chemistry/$(CHMDIR)/second_Main.f90 second_Model.o second_Initialize.o
	 $(F90) $(CHEM_OPTS) -c $< -o $@

# ACDC, NH3
$(OBJDIR)/%.o: $(SRCDIR)/ACDC/ACDC_module_ions_2018_08_31/%.f90
	$(F90) $(ACDC_OPTS) -c $< -o $@

# ACDC, DMA
$(OBJDIR)/%.o: $(SRCDIR)/ACDC/ACDC_module_2016_09_23/%.f90
	$(F90) $(ACDC_OPTS) -c $< -o $@


# Common solvers for ACDC
$(OBJDIR)/solution_settings.o: ACDC/ACDC_module_ions_2018_08_31/solvers/solution_settings.f90
	 $(F90) $(ACDC_OPTS) -c $< -o $@
#
$(OBJDIR)/vode.o: ACDC/ACDC_module_ions_2018_08_31/solvers/vode.f
	$(F90) -std=legacy -O3 -c $< -o $@

$(OBJDIR)/vodea.o: ACDC/ACDC_module_ions_2018_08_31/solvers/vodea.f
	$(F90) -std=legacy -O3 -c $< -o $@

$(OBJDIR)/chemistry.o: chemistry.f90 second_Parameters.o second_Main.o
	$(F90) $(BOX_OPTS) -c $< -o $@

# Actual model files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(F90) $(BOX_OPTS) -c $< -o $@ $(NETLIBS)

BOX_MODS = $(BOX_OBJECTS:.o=.mod)

CHEM_MODS = (second_precision.mod second_monitor.mod second_parameters.mod second_initialize.mod second_util.mod second_jacobiansp.mod \\
                second_linearalgebra.mod second_jacobian.mod second_global.mod second_rates.mod second_integrator.mod second_function.mod \\
                second_model.mod second_main.mod)

ACDC_MODS = $(ACDC_OBJECTS:.o=.mod)
ACDC_D_MODS = $(ACDC_D_OBJECTS:.o=.mod)

# This entry allows you to type 'make clean' to get rid of all object and module files
# With 'clean', don't remove chemistry object files, since it takes very long (30 min.) to compile them,
# and there usually is no need to recompile them

clean:
	-@rm $(BOX_OBJECTS) $(BOX_MODS) 2>/dev/null || true
	-@cd $(OBJDIR) ; rm arcabox.o   2>/dev/null || true
	-@rm arcabox.exe                2>/dev/null || true

cleaneverything:
	-@cd $(OBJDIR) ; rm -r *        2>/dev/null || true
	-@rm arcabox.exe                2>/dev/null || true

clean_current_chemistry:
	-@rm $(BOX_OBJECTS) $(BOX_MODS) 2>/dev/null || true
	-@cd $(OBJDIR)/$(CHMDIR) ; rm *.mod *.o    2>/dev/null || true
	-@cd $(OBJDIR) ; rm arcabox.o    					 2>/dev/null || true
	-@rm arcabox.exe                					 2>/dev/null || true


""")
f.close()

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
    f.write('%s'%le)
    if csc:
        f.write('module load python-data/3.7.6-1%s'%le)
        f.write('module load gcc/9.1.0%s'%le)
        f.write('module load netcdf/4.7.0%s'%le)
        f.write('module load netcdf-fortran/4.4.4%s'%le)
    else:
        print('ARCA user interface can be used on CSC with NoMachine. To get correct settings,')
        print('run this script with --csc flag or edit run_arca.sh and uncomment the csc options.')
        f.write('# If you are using the user interface on CSC, uncomment next 4 lines%s'%le)
        f.write('# module load python-data/3.7.6-1%s'%le)
        f.write('# module load gcc/9.1.0%s'%le)
        f.write('# module load netcdf/4.7.0%s'%le)
        f.write('# module load netcdf-fortran/4.4.4%s'%le)

    f.write('%s'%le)
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
