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

dt = 10.0 

true  ='.True.'
false ='.False.'
### dimensions of the file

file_check = np.loadtxt(path_cases_dir+path_case_name+'/'+Fname)
(rr,col)=file_check.shape

## write to this file
file = open("/home/local/carltonx/Box-model/supermodel-phase-1/model_init_test" , "w+")

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
file.write("dt   =  %f\n" %dt  )
file.write("/ \n\n")

file.write(" ####### file size ##### \n\n")

file.write("&NML_shape \n ")
file.write("rows   =  %d \n " %rr)
file.write("cols   =  %d\n" %col)
file.write("/ \n\n")

file.close()
