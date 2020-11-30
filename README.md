You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too. ON WINDOWS, THE BEST WAY INSTALL THE NECESSARY PACKAGES IS WITH CYGWIN.
See this tutorial video: https://www.helsinki.fi/fi/unitube/video/7419bbe3-3fc8-493a-b076-7307d2e8191c

To use the graphical user interface, you need Python 3. We strongly encourage to use the native Python on Windows
instead of Linux Subsystem (although this may work also).

After Python 3 is available, there is a client to set up the model. Go to directory install/, start the terminal and run
setup.py (for example, if you call Python 3 with python3, the command is "python3 setup.py"). The setup will ask you to
install the necessary Python Modules, if they are not available. If you want to use some other installation procedure for
them, just answer "n" for the prompt.

Next you will be asked if the Fortran model should be compiled. If this is succesful*, a script file is written in the
root folder, and the GUI can be started from there by calling "sh run_arca.sh"


The necessary Python modules are:

- numpy

- scipy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)

- matplotlib (also optional, but recommended)




Learn more of using ARCA from this video: https://www.helsinki.fi/fi/unitube/video/b2c6775e-ad5f-4cee-9fd5-e502af2fb256

Note that the program is constantly developed and the look and functionality may differ from this video.


* If the compiling failed, it might because the netCDF4 libs were not properly configured in the makefile. These would
vary from one system to another, and therefore we cannot guarantee that the makefile currently is able to find them. In
this case open the makefile, and comment the line which starts "NETLIBS" and uncomment the next line. Also if you are
using some other Fortran compiler than gfortran, you need to define this in the beginning of the makefile (F90 = gfortran)
