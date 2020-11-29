To use the model, please refer to Documentation.

You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too.

For convenience, there is a client to set up and run the model, this is written in Python 3.6. Necessary libraries are:

- numpy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)


You can (probably) install these with command: <python 3> -m pip install -U <module>, but this might vary from one
Python installation to another. Also using Python VIRTUALENV is advisable.

TO USE THE CLIENT:
After the python modules are working, to setup the client, open therun.sh in text editor and change the path after "cd"
in line 2 to the directory where this readme-file and makefile are located.

For example, if this file is located in: "~/models/thebox", change the line to

cd ~/models/thebox

Run the client by executing thebox.sh
