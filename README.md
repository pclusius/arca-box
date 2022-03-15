Atmospherically Relevant Chemistry and Aerosol box model
Copyright (C) 2021  Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
Contact information arca@helsinki.fi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

----------------------------------------------------------------------


How to Install ARCA and the necessary software for it
-----------------------------------------------------
You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too. The user interface is written in Python 3.

How to install on Linux Mint / Ubuntu:
---------------------------------------

https://www.helsinki.fi/fi/unitube/video/80d7a662-f60d-4cd1-b648-31dd0cbe9bd2


How to install on Puhti/CSC:
----------------------------
1) Go to directory install/
2) run install_on_csc.sh (e.g. "sh install_on_csc.sh")
3) follow the instructions. If the compiling step does not work for some reason, you can try running the install script
   again but omit the compiling. Then navigate to ARCA root dir and try compiling there again (with command "make").
   In the "makefile", you might want to change the code optimization to O2 or even O3 (default is O1).


How to install on Windows 10:
-----------------------------

The best way install the necessary packages is with Cygwin. See this tutorial video:
https://www.helsinki.fi/fi/unitube/video/7419bbe3-3fc8-493a-b076-7307d2e8191c

To use the graphical user interface, you need Python 3. We strongly encourage to use the native Python on Windows
instead of Linux Subsystem (although this may work also).


For all systems (except Puhti/CSC), when Fortran and Python 3 is available:
---------------------------------------------------------------------------

After Python 3 is available, there is a script to set up the model. Go to directory install/, start the terminal and run
setup.py (for example, if you call Python 3 with python3, the command is "python3 setup.py"). The setup will ask you to
install the necessary Python Modules, if they are not available. If you want to use some other installation procedure for
them, or you know you already have the necessary modules, just answer "n" for the prompt.

Next you will be asked if the Fortran model should be compiled. If this is succesful*, a starter script file is written
in the root folder, and the GUI can be started from there by calling "sh run_arca.sh"


\* If the compiling failed, it might because the netCDF4 libs were not properly configured in the makefile. These would
vary from one system to another, and therefore we cannot guarantee that the makefile currently is able to find them. In
this case open the makefile, and comment the line which starts "NETLIBS" and uncomment the next line. Also if you are
using some other Fortran compiler than gfortran, you need to define this in the beginning of the makefile (F90 = gfortran)


The necessary Python modules are:
---------------------------------

- numpy

- scipy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)

- matplotlib (also optional, but recommended)

The setup.py script tries to install these using command [python3] -m pip install \--user [module], where [python3] is
the call for the user's Python 3.




How to use ARCA box?
--------------------
Learn more of using ARCA from this video: https://www.helsinki.fi/fi/unitube/video/b2c6775e-ad5f-4cee-9fd5-e502af2fb256

Note that the program is constantly developed and the look and functionality may differ from this video.



Installation troubleshooting
----------------------------
If you get errors when installing the Python modules, try upgrading pip:

python3 -m pip install --upgrade pip

and then run the setup.py again.


Version history
---------------

1.2.0

What's new:
- 5 ACDC systems that are dynamically used so that their monomer names and system size can be changed without hard
  coding, needs recompiling though
- ACDC systems have proper recording of the inut data files and system names
- Added new section to the GUI->"Cluster formation"-tab, where ACDC components can be linked with model components
- Updated ACDC plugin
- Hyde3 chemistry scheme replaced Hyde2. Inorganic precursors are fixed (their concentration time derivatives are 0).
  This steady state approach is already used in some OVOCs in the scheme. The user is advised to check the compounds in
  listed unde DEFFIX in file second.def and use the appropriate approach for their work.
- Drag and drop INITFILE loading
- Losses calculated with DVODE
- loss parametrization clarified
- Updated spectral data
- It is possible to send in direct actinic flux data
- (multi)modal distribution shows mass and total area
- ARCA paper submission version

Fixes:
- Fixed SMEAR II short wave radiation spectrum
- Windows installer now prefers pyqtgraph 0.12.0
- Minor bug fixes in GUI and Fortran model


1.1.2
Hotfix:
- fixed error in calculation of Koehler factor for particles of size 1 nm (in aerosol_dynamics.f90)


1.1.1
Fixes:
- Fixed empty array problem in second_reactivity.90 and add_reactivity.py (producing error with some fortran compilers).
- Fixed makefile (version control gave problems on some Win systems)


1.1.0

ARCA is now licenced under GNU GPL


1.0.6

Fixes:
- Kelvin term was approximated using Taylor series. This is a bad approximation with very small particles.
  Now Kelvin term is calculated using exponential form. To use the old Kelvin term (used by many other models and the
  ACP scheme), set custom option (in NML_CUSTOM or in the GUI Run ARCA -> Custom model options) Kelvin_taylor = .true.

(Further additions in current version, will be put to next release)
  What's new:
  Fixes:


1.0.5

What's new:
- makefile supplemented with Puhti configurations
- Output directories are now created automatically

Fixes:
- makefile improved; compiling chemistry module for the first time should not lead to error
- On mac, the setup.py creates run_arca.command instead of run_arca.sh; is also be double-clickable
- gui: replaced some deprecated Pyqtgraph methods giving errors on pyqtgraph 0.12
- gui: fixed duplicating legend in mass plot


1.0.4

What's new:
- added GPL licence to KPP and ACDC directories.
- added coag_sink to Particles.nc
- added CS_calc to General.nc

Fixes:
- Added new Hyde chemistry to accommodate reactivity calculations
- fixed color issues in mass plotting
- fix for how losses file was interpolated


1.0.3

What's new:
- Create second_reactivity.f90
Fixes:
- minor fixes in the mass plotting
