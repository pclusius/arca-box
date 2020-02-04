To use the model, please refer to Documentation.

You will need C and Fortran compilers (gcc, gfortran, iFortran and/or similar) to compile the model. Model relies on
NetCDF4 to save files, so you'll need that too.

For convenience, there is a client (run it by executing superbox_client.sh) to set up and run the model, this is written in Python 3.6. Necessary libraries are:

- numpy

- PyQt5

- pyqtgraph

- netcdf4-python (optional, but highly recommended)

You can probably install these with command: <python 3> -m pip install -U <module>, but this might vary from one
Python installation to another. Also using Python VIRTUALENV is advisable.
