To use the model, please refer to Documentation.

You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too.

For convenience, there is a client to set up and run the model, this is written in Python 3.6. Necessary libraries are:

- numpy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)

To run the client, open therun.sh in text editor and change the path after "cd" in line 2 to the directory where this
readme-file and makefile are located:

for example, if this file is located in: ~/05-APCAD/supermodel-phase-1, change the line to

cd ~/05-APCAD/supermodel-phase-1

Run the client by executing thebox.sh

You can (probably) install these with command: <python 3> -m pip install -U <module>, but this might vary from one
Python installation to another. Also using Python VIRTUALENV is advisable.
