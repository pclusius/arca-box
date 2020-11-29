You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too. ON WINDOWS, THE BEST WAY TO THIS IS WITH CYGWIN.

For convenience, there is a client to set up and run the model, this is written in Python 3.6. Necessary libraries are:

- numpy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)


You can (probably) install these with command: <python 3> -m pip install -U <module>, but this might vary from one
Python installation to another. Also using Python VIRTUALENV is advisable.

TO USE THE CLIENT:
After the python modules are working, to setup the client, run setup.py and see that the model compiles. 
The script asks for the command used to run python, enter it. After the script is finished, start the model 
by typing the command: 

sh run_arca.sh
