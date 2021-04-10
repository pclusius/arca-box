#!/bin/bash

module load gcc/9.1.0
module load netcdf/4.7.0
module load netcdf-fortran/4.4.4

python3 setup.py --csc
