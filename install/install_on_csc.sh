#!/bin/bash


cluster=$1

case "$cluster" in
    puhti)
        module load python-data
        module load gcc
        module load netcdf-c
        module load netcdf-fortran
        ;;
    mahti)
        module load python-data
        module load gcc
        module load netcdf-c
        module load netcdf-fortran
        ;;
    lumi)
        module purge
        module load cray-python
        module load PrgEnv-gnu
        module load cray-hdf5
        module load cray-netcdf
        module load craype-x86-rome  # target architecture of login node
        ;;
    *)
        echo "Usage: $0 {puhti|mahti|lumi}"
        exit 1
        ;;
esac

python3 setup.py --csc "$cluster"