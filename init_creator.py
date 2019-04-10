#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C#for creating a .init file
"""
import os
import numpy as np
import sys
import pylab as pl 


path_dir           = '/home/local/carltonx/Box-model/supermodel-phase-1/'
path_cases_dir     = '/home/local/carltonx/Box-model/supermodel-phase-1/Cases/'
path_case_name     = 'Test_case_01'
Fname              = 'Api_100.txt'
Dfile              = 'dmps.txt'
Temp_file          = 'SMEAR_TEMP_168.dat'
path_pec           = '/home/local/carltonx/Box-model/PC180330/'
Voc_file           = 'VOC_0042.dat'

dt       = 10.0  ## seconds
sim_time = 3600  ## seconds
ctime    = 0.0   ## current time. 0 at beginning

true  ='.True.'
false ='.False.'
### dimensions of the file

file_check = np.loadtxt(path_cases_dir+path_case_name+'/'+Fname)
dmps_check = np.loadtxt(path_cases_dir+path_case_name+'/'+Dfile)
temp_check = np.loadtxt(path_pec + Temp_file)
voc_check  = np.loadtxt(path_pec + Voc_file )

(rr,col)      =file_check.shape
(rr_d, col_d) =dmps_check.shape
(rr_t, col_t) =temp_check.shape
(rr_v, col_v) =voc_check.shape

## write to this file
file = open("/home/local/carltonx/Box-model/supermodel-phase-1/model_init" , "w+")

##
file.write("####### model initialization setup ############\n\n")

file.write("&NML_path \n ")
file.write("Work_Dir    =  '%s' \n " %path_dir  )
file.write("Case_Dir    =  '%s' \n " %path_cases_dir  )
file.write("Case_name   =  '%s' \n " %path_case_name  )
file.write("Input_file  =  '%s'\n" %Fname  )
file.write("/ \n\n")
 


file.write("####### Flags setup True or false ####### \n\n")

file.write("&NML_Flag \n ")
file.write("Aerosol_flag   =   %s \n " %true  )
file.write("Chemistry_flag =   %s \n " %true  )
file.write("Particle_flag  =   %s \n " %false  )
file.write("Extra_data     =   %s\n" %false  )
file.write("/ \n\n")

file.write("####### Time flag ######## \n\n")

file.write("&NML_Time \n ")
file.write("dt         =  %f \n " %dt  )
file.write("sim_time   =  %f \n " %sim_time )
file.write("ctime      =  %f\n" %ctime )
file.write("/ \n\n")

file.write(" ####### file size ##### \n\n")

file.write("&NML_shape \n ")
file.write("rows   =  %d \n " %rr)
file.write("cols   =  %d\n" %col)
file.write("/ \n\n")

file.write(" ####### dmps file read ##### \n\n")

file.write("&NML_DMPS \n ")

file.write("DMPS_file         =  '%s' \n " %Dfile)
file.write("dmps_file_check   =   %s \n "  %true)
file.write("dmps_rows         =   %d \n " %rr_d)
file.write("dmps_cols         =   %d\n" %col_d)
file.write("/ \n\n")

file.write(" ####### TEMP file read ##### \n\n")

file.write("&NML_TEMP \n ")
file.write("TEMP_path         =  '%s' \n " %path_pec)
file.write("TEMP_file         =  '%s' \n " %Temp_file)
file.write("temp_file_check   =   %s \n "  %true)
file.write("temp_rows         =   %d \n " %rr_t)
file.write("temp_cols         =   %d\n" %col_t)
file.write("/ \n\n")


file.write(" ####### VOC file read ##### \n\n")

file.write("&NML_VOC \n ")
file.write("VOC_file          =  '%s' \n " %Voc_file)
file.write("voc_file_check   =   %s \n "  %true)
file.write("voc_rows         =   %d \n " %rr_v)
file.write("voc_cols         =   %d\n" %col_v)
file.write("/ \n\n")

file.close()
