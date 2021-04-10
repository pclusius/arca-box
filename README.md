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
