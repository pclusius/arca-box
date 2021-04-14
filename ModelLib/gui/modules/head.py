# -*- coding: utf-8 -*-
"""
=============================================================================
Copyright (C) 2021  Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
Contact information arca@helsinki.fi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================
"""

## Some constants --------------------------------------------
# widths of the columns in "Input variables" tab
column_widths = [120,90,70,70,70,90,50,3]

# available units for variables, used to fill the tables and graphs with appropriate units
units = {
'TEMPK': ['K','C'],
'TIME_IN_HRS': ['hrs'],
'TIME_IN_SEC': ['s'],
'PRESSURE': ['Pa','hPa','mbar','kPa','bar','atm'],
'REL_HUMIDITY': ['%'],
'CONDENS_SINK':['1/s'],
'CON_SIN_NITR':['1/s'],
'SW_RADIATION':['W/m2'],
'ION_PROD_RATE':['ip/cm3 s'],
'NUC_RATE_IN':['1/cm3 s'],
'REST':['#/cm3','ppm','ppb','ppt','ppq']
}
# Name of the executable -------------------------------------------
exe_name = 'arcabox.exe'
# Path to variable names -------------------------------------------
path_to_names = 'ModelLib/required/NAMES.dat'
# GUI root
gui_path = 'ModelLib/gui/'
# GUI defaults are saved into this file. If it doesn't exist, it gets created in first start
defaults_file_path = gui_path+'conf/defaults.init'
minimal_settings_path = gui_path+'conf/minimal.init'
monitorfont_pickle = 'conf/monitorfont.pickle'
globalfont_pickle = 'conf/globalfont.pickle'

# This path will be added to Common out if no other option is given
default_inout = 'INOUT'
# This case name will be used as case name if no other name is given
default_case = 'DEFAULTCASE'.upper()
# This run name will be used as run name if no other name is given
default_run  = 'DEFAULTRUN'.upper()
# name and location of the temporary settings file used for test runs
tempfile = gui_path+'tmp/GUI_INIT.tmp'
# initial maximum for function creator tab sliders
slMxs = (200,190,220,100,200)
# 10 colors for plots
colors = [(120,0,0),(180,0,0),(220,0,0),(255,10,0),(255,85,0),
(255,85,140),(255,85,255),(180,0,255),(110,0,255),(0,0,255)]
# BG colour for ENV vars
env_no = (215,238,244)
env_yes = (128,179,255)
# BG colour for ORG vars
org_no = (243,240,239)
org_yes = (172,147,147)

# icon
modellogo = gui_path+"/icons/ArcaLogo.png"
boxicon = gui_path+"/icons/thebox_ico.png"
CurrentVersion = "ARCA Box Model 0.9"
# Some messages
netcdfMissinnMes = ('Please note:',
'To open NetCDF-files you need netCDF4 for Python.\nYou can istall it with pip, package manager (or perhaps: python3 -m pip install --user netCDF4.')


# files that can be modified with the Editor
nucl_homs = "ModelLib/required/nucl_homs.txt"
custom_functions = "src/custom_functions.f90"
AmmSystemFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/input_ANnarrow_neutral_neg_pos.inp"
Amm_EnergyFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/HS298.15K_426clusters2016Apr25.txt"
Amm_DipoleFile = "src/ACDC/ACDC_module_ions_2018_08_31/Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt"
DMASystemFile = "src/ACDC/ACDC_module_2016_09_23/Perl_input/input_AD.inp"
DMA_EnergyFile = "src/ACDC/ACDC_module_2016_09_23/Perl_input/dH_dS.txt"
SCREENPRINT_NML = "ModelLib/required/NML_SCREENPRINTS.def"
# Create chemistry script location
ccloc = 'ModelLib/gui/chemistry_package_PZ'

# Guess initially the current python version and check the calling script for precise version
currentPythonVer = 'python'
try:
    with open('run_arca.sh') as f:
        for line in f:
            if gui_path in line:
                currentPythonVer = line.split()[0]
except:
    pass
