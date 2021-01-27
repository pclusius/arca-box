#!/usr/bin/env python3

#==============================================================================#
#
# Created by Zhou Putian, University of Helsinki
#
# Combination part refers to Sampo Smolander's code
#
#==============================================================================#


#==============================================================================#
#
# Header
#
#==============================================================================#
import os
import sys

import logging
import datetime

import re
import argparse


#==============================================================================#
#
# Helper functions
#
#==============================================================================#

# Find first index containing substring
def index_containing_substring(the_list, substring, case_sensative=True):
  for i, s in enumerate(the_list):
    if case_sensative:
      mat = re.search(substring, s)
    else:
      mat = re.search(substring, s, re.I)

    if mat:
      return i, mat.start(), mat.end()

  return -1, -1, -1

def update_defvar_list(file_lines, species_list, defvar_print_list):
  # Find '#DEFVAR'
  line_number_defvar, i1, i2 = index_containing_substring(file_lines, r'^\s*#DEFVAR', False)

  # Start only when there is a #DEFVAR block
  if line_number_defvar >= 0:
    # Start from the next line
    for i, s in enumerate(file_lines[line_number_defvar+1:]):
      # Meet next block '#XXXX'
      if re.search(r'^\s*#', s):
        return
      # Still in the block '#DEFVAR'
      else:
        # Search for, e.g., 'OH = IGNORE ;'
        # \S: not white space, opposite to \s
        mat = re.search(r'^\s*(\S+)\s*=\s*\S+.*;', s)
        # Found species
        if mat:
          # Not duplicated, case insensitive
          if mat.group(1).upper() not in (tmp.upper() for tmp in species_list):
            species_list.append(mat.group(1))
            defvar_print_list.append(s)
        # Found other statement
        else:
          defvar_print_list.append(s)
        
  return


def update_deffix_list(file_lines, species_list, deffix_print_list):
  # Find '#DEFFIX'
  line_number_deffix, i1, i2 = index_containing_substring(file_lines, r'^\s*#DEFFIX', False)

  # Start only when there is a #DEFFIX block
  if line_number_deffix >= 0:
    # Start from the next line
    for i, s in enumerate(file_lines[line_number_deffix+1:]):
      # Meet next block '#XXXX'
      if re.search(r'^\s*#', s):
        return
      # Still in the block '#DEFFIX'
      else:
        # Search for, e.g., 'OH = IGNORE ;'
        # \S: not white space, opposite to \s
        mat = re.search(r'^\s*(\S+)\s*=\s*\S+.*;', s)
        # Found species
        if mat:
          # Not duplicated, case insensitive
          if mat.group(1).upper() not in (tmp.upper() for tmp in species_list):
            species_list.append(mat.group(1))
            deffix_print_list.append(s)
        # Found other statement
        else:
          deffix_print_list.append(s)
        
  return


def update_ro2_list(file_lines, ro2_list):
  # Find the starting line of "RO2 = "
  line_number_ro2, i1, i2 = index_containing_substring(file_lines, r'^\s*RO2\s*=', False)

  # Start only when there is a RO2 definition
  if line_number_ro2 >= 0:
    # Starting to count RO2 species
    for i, s in enumerate(file_lines[line_number_ro2+1:]):
      # Skip empty lines
      if re.search(r'^\s*$', s):
        continue
      # RO2 definition ending
      elif not re.search(r'C\(ind_', s, re.I):
        return
      # Still inside the RO2 definition
      else:
        # Split by '+'
        c_list = s.split('+')
        for c in c_list:
          # Search for RO2 species
          mat = re.search(r'C\(\s*ind_(\S+)\s*\)', c, re.I)
          # Found species
          if mat:
            # Not duplicated
            if mat.group(1) not in ro2_list:
              ro2_list.append(mat.group(1))

  return 


def update_equation_list(file_lines, equation_list, equation_norate_list, equation_print_list, equation_duplicate_print_list):
  # Find '#EQUATIONS'
  line_number_equations, i1, i2 = index_containing_substring(file_lines, r'^\s*#EQUATIONS', False)

  # Start only when there is a #EQUATIONS block
  if line_number_equations >= 0:
    # Start from the next line
    for i, s in enumerate(file_lines[line_number_equations+1:]):
      # Meet next block '#XXXX'
      if re.search(r'^\s*#', s):
        return
      # Still in the block '#EQUATIONS'
      else:
        # Search for, e.g., "{10.}    NO + NO3 = NO2 + NO2 :   1.8D-11*EXP(110/TEMP)   ;"
        if re.search('\{', s):  # with {tag} in the beginning
          mat = re.search(r'^\s*{.*}\s*(\S+.*)=(\s*\S+.*):(\s*\S+.*);', s)
        else:  # without {tag}
          mat = re.search(r'^\s*(\S+.*)=(\s*\S+.*):(\s*\S+.*);', s)
        # Found equations
        if mat:
          # Reactants, products and rate coefficients without spaces
          reactants = [x.strip() for x in mat.group(1).split('+')]
          products = [x.strip() for x in mat.group(2).split('+')]
          ratecoef = mat.group(3).strip()

          # Whole equation with sorted reactants and products
          equation = []
          equation.extend(sorted(reactants))
          equation.append('=')
          equation.extend(sorted(products))
          equation.append(':')
          equation.append(ratecoef)

          # print(equation)

          # Equation with sorted reactants and products and without rate coefficient
          equation_norate = []
          equation_norate.extend(sorted(reactants))
          equation_norate.append('=')
          equation_norate.extend(sorted(products))

          # print(equation_norate)

          # Not duplicated equation without considering rate
          if equation_norate not in equation_norate_list:
            equation_list.append(equation)
            equation_norate_list.append(equation_norate)
            equation_print_list.append(s)
            # print(equation_norate)
          # Duplicated equation without rate
          else:
            # Only the rates are different, then save the line to equation_duplicate_print_list
            # Still print out to final modified file, but keep the information in a log file
            if equation not in equation_list:
              equation_list.append(equation)
              equation_norate_list.append(equation_norate)
              equation_print_list.append(s)
              equation_duplicate_print_list.append(s)
        # Found other statement, e.g., comments starting with //
        else:
          equation_print_list.append(s)

  return


if __name__ == '__main__':
  #==============================================================================#
  #
  # Input variables
  #
  # 1. Set the name of originally MCM generated KPP file: mcm_kpp_file_name.
  # 2. mcm_molar_mass_file_name is not used now.
  # 3. Set the name of final KPP def file modified_mcm_kpp_file_name,
  #    usually it is 'second.def'.
  # 4. Set the log file name log_file_name to save the running information for
  #    debugging.
  # 5. Set the inline def file name inline_file_name, usually it is 'inline.def'
  #    or 'File_1_inline.def'.
  # 6. If you want to combine #DEFVAR and #EQUATIONS from other files, you can
  #    add their names to the list include_file_list.
  # 7. Set is_ciso3 to True if you want to update reaction rates of CI -> SO3
  #    from 7.00D-14 to 7.00D-13.
  # 8. Run the code as:
  #    $ python create_chemistry.py
  #
  # Future work
  # 1. Deal with multiple declaration for species.
  # 2. Show which files the duplicated equations are from.
  # 3. Get arguments as input for this python script.
  # 4. NA --> HNO3
  # 5. Ask if the files are included.
  # 6. Ask if start kpp now.
  # 7. Calculate vapor pressure.
  #
  #==============================================================================#
  
  parser = argparse.ArgumentParser(description='Creating new chemistry files.', add_help=True)
  subparsers = parser.add_subparsers(help='Commands')
  
  kppdef_parser = subparsers.add_parser('kppdef', help='Create modified second.def after combination.')
  kppdef_parser.add_argument('mcm_kpp_file_name', default='mcm_subset.kpp', action='store', \
      help='Absolute path of the input KPP file generated by MCM.')
  kppdef_parser.add_argument('-o', dest='modified_mcm_kpp_file_name', default='second.def', action='store', \
      help='Absolute path of the modified MCM KPP file.')
  kppdef_parser.add_argument('--log', '-l', dest='log_file_name', default='create_chemistry.log', action='store', \
      help='Log file name.')
  kppdef_parser.add_argument('--include-files', '-f', \
      dest='include_file_list', \
      nargs='+', \
      default=['File_2_PRAM_v21.txt', 'File_3_terpenes_not_in_mcm.txt', 'File_4_emission.txt'], \
      action='store', \
      help='A list of included files. Specify like: -f a.txt b.txt ...')

  args = parser.parse_args()

  # parser.add_argument('integers', metavar='N', type=int, nargs='+',
  #                     help='an integer for the accumulator')
  # parser.add_argument('--sum', dest='accumulate', action='store_const',
  #                     const=sum, default=max,
  #                     help='sum the integers (default: find the max)')
  # 
  # args = parser.parse_args()
  # print(args.accumulate(args.integers))
  # 
  # sys.exit()
  
  # Absolute path of the input KPP file generated by MCM
  # mcm_kpp_file_name = './mcm_subset.kpp'
  
  # Absolute path of the molar mass file generated by MCM
  # Columns after header: name, SMILES, InChI, molar mass (g mol-1)
  # Inorganic species are not included
  mcm_molar_mass_file_name = './mcm_subset_mass.txt'
  
  # Absolute path of the modified MCM KPP file
  # modified_mcm_kpp_file_name = './second.def'
  
  # Log file name
  # log_file_name = 'create_chemistry.log'
  
  # The file name of inline.def
  inline_file_name = 'File_1_inline.def'
  
  # Included file list
  include_file_list = [s.strip() for s in args.include_file_list]  # strip the white spaces
  # include_file_string_list = ['#INCLUDE ' + x + '\n' for x in include_file_list]
  print(args.include_file_list)
  print(include_file_list)
  
  # H2SO4 = dummy : RES1
  string_h2so4_dummy = 'H2SO4 = dummy : RES1 ;\n'
  
  # If Criegee intermediate (CI) radicals with e.g. SO3 will be considered
  is_ciso3 = True
  
  
  #==============================================================================#
  #
  # Initiate log file
  #
  #==============================================================================#
  logging.basicConfig( \
    level = logging.DEBUG, \
    filename = args.log_file_name, \
    filemode = 'w', \
    format = '%(message)s' \
    )
  
  logging.info( 'Created at ' + datetime.datetime.now().strftime("%H:%M:%S, %B %d, %Y") )
  
  
  #==============================================================================#
  #
  # Fix text bugs in mcm_kpp_file_name
  #
  #==============================================================================#
  
  # Read the whole file and save it to a string list
  with open(args.mcm_kpp_file_name, 'r') as f:
    file_lines = f.readlines()
  
  
  
  #---------- Add the complex reaction rates ----------#
  # Add the content of File_1_inline_def to the root
  # def file.
  #----------------------------------------------------#
  
  if os.path.isfile(inline_file_name):
    # Get the file content of inline def file
    with open(inline_file_name, 'r') as f:
      inline_file_lines = f.readlines()
    
    # Insert the content after '#INCLUDE atoms'
    line_number_include_atoms, i1, i2 = index_containing_substring(file_lines, r'#include atoms', False)
    file_lines[line_number_include_atoms+1:line_number_include_atoms+1] = ['\n'] + inline_file_lines + ['\n']
    
    # Write to log file
    logging.info('')
    logging.info('Added complex reaction rates from {0}.'.format(inline_file_name))
  else:
    # Write to log file
    logging.info('')
    logging.info('{0} does not exist.'.format(inline_file_name))
  
  
  #---------- Include files ----------#
  # Write the '#INCLUDE files' after '#INCLUDE atoms'.
  # No needed now, since all the #DEFVAR and #EQUATIONS are already included
  # in the code below.
  #-----------------------------------#
  
  # file_lines[line_number_include_atoms+1:line_number_include_atoms+1] = ['\n'] + include_file_string_list + ['\n']
  
  
  #---------- Fix bugs ----------#
  
  # = IGNORE --> dummy = IGNORE
  line_number_empty_ignore, i1, i2 = index_containing_substring(file_lines, r'^\s*=\s*IGNORE', False)
  if line_number_empty_ignore >= 0:
    file_lines[line_number_empty_ignore] = 'dummy = IGNORE ;\n'
  
    # Write to log file
    logging.info('')
    logging.info('Fix the bug: "= IGNORE --> dummy = IGNORE".')
  
  # Delete 'USE constants'
  line_number_use_constants, i1, i2 = index_containing_substring(file_lines, r'\s*USE\s*constants', False)
  if line_number_use_constants >= 0:
    file_lines[line_number_use_constants] = ''
  
    # Write to log file
    logging.info('')
    logging.info('Deleted the line "USE constants".')
  
  # Delete 'Call mcm_constants(..)'
  line_number_call_constants, i1, i2 = index_containing_substring(file_lines, r'\s*CALL\s*mcm_constants', False)
  if line_number_call_constants >= 0:
    file_lines[line_number_call_constants] = ''
  
    # Write to log file
    logging.info('')
    logging.info('Deleted the line "CALL mcm_constants(time, temp, M, N2, O2, RO2, H2O)".')
  
  # "O + O3 = xxx :" --> "O + O3 = dummy :"
  # Notice: double check is used here to avoid changing, e.g., CH3O + O3 = xxx.
  line_number_o_o3, i1, i2 = index_containing_substring(file_lines, r'O\s*\+\s*O3\s*=\s*:', False)
  if line_number_o_o3 >= 0:
    tmp_string = file_lines[line_number_o_o3]
    file_lines[line_number_o_o3] = tmp_string[:i1] + 'O + O3 = dummy :' + tmp_string[i2:]
  
    # Write to log file
    logging.info('')
    logging.info('Fix the bug: "O + O3 = --> O + O3 = dummy".')
  
  # "OH + HO2 = :" --> "OH + HO2 = dummy :"
  line_number_oh_ho2, i1, i2 = index_containing_substring(file_lines, r'OH\s*\+\s*HO2\s*=\s*:', False)
  if line_number_oh_ho2 >= 0:
    tmp_string = file_lines[line_number_oh_ho2]
    file_lines[line_number_oh_ho2] = tmp_string[:i1] + 'OH + HO2 = dummy :' + tmp_string[i2:]
  
    # Write to log file
    logging.info('')
    logging.info('Fix the bug: "OH + HO2 = --> OH + HO2 = dummy".')
  
  
  #---------- Stuff related to H2SO4 ----------#
  
  # "SA" --> "H2SO4"
  count_SA = 0
  for i, s in enumerate(file_lines):
    if re.search(r'\bSA\b', s, re.I):
      count_SA += 1
      file_lines[i] = re.sub(r'\bSA\b', r'H2SO4', s)
  
  # Write to log file
  logging.info('')
  logging.info('{0} SA are changed to H2SO4.'.format(count_SA))
  
  
  # Add "H2SO4 = dummy : RES1"
  file_lines.append(string_h2so4_dummy)
  
  # Write to log file
  logging.info('')
  logging.info('Add condensation sink RES1 for H2SO4: "{0}"'.format(re.sub(r'\n$', '', string_h2so4_dummy)))
  
  
  #---------- Comment out the unused reaction rate coefficients calculation ----------#
  # KNO = KRO2NO*NO
  # KHO2 = KRO2HO2*HO2*0.706
  # KRO2 = 1.26D-12*O2
  # KNO3 = KRO2NO3*NO3
  # KTR = KNO + KHO2 + KRO2 + KNO3
  # K16ISOM = (KTR*5.18D-04*EXP(1308/TEMP)) +(2.76D07*EXP(-6759/TEMP))
  #-----------------------------------------------------------------------------------#
  
  # Write a newline in log file
  # logging.info('')
  
  for str_k in ['KNO', 'KHO2', 'KRO2', 'KNO3', 'KTR', 'K16ISOM']:
    line_number_k, i1, i2 = index_containing_substring(file_lines, str_k + r'\s*=', False)
    if line_number_k >= 0:
      file_lines[line_number_k] = '! ' + file_lines[line_number_k]
  
      # Write to log
      logging.info('Comment out {0}'.format(str_k))
  
  
  #==============================================================================#
  #
  # Combine different include files to one file
  #
  # 1. Find all the include files written in the root def file
  # 2. Combine chemical species defined for #DEFVAR
  # 3. Combine chemical species defined for #DEFFIX
  # 4. Combine RO2 species
  # 5. Combine chemical reactions for #EQUATIONS
  #
  #==============================================================================#
  
  #---------- 1. Find all the include files in the root def file but do not count atoms ----------#
  
  # Write to log file
  logging.info('')
  
  # File all the files from '#INCLUDE' lines
  for i, s in enumerate(file_lines):
    match_file_name = re.search(r'^\s*#INCLUDE\s*(\S+)\s*$', s, re.I)
    if match_file_name:
      tmp_file_name = match_file_name.group(1)
      if tmp_file_name not in include_file_list:
        include_file_list.append(match_file_name.group(1))

  # Do not count 'atoms'
  if 'atoms' in include_file_list:
    include_file_list.remove('atoms')
  
  # Check if include files exist
  tmp_list = list(include_file_list)  # avoiding errors in the list iteration
  for f in include_file_list:
    if not os.path.isfile(f):
      # Remove from the list if not existing
      tmp_list.remove(f)
  
      # Write the non-existing file name to log file
      logging.info('{0} does not exist.'.format(f))

  include_file_list = list(tmp_list)
  
  # Write included files to log file
  logging.info('The included files besides the root def file (the inline file {0} has already been included):'.format(inline_file_name))
  for f in include_file_list:
    logging.info(f)


  #---------- 1.5. Add #INLINE parts (e.g., F90_GLOBAL, F90_RCONST) ----------#

  # print(include_file_list)
  
  
  #---------- 2. Combine chemical species defined for #DEFVAR ----------#
  
  # Initiate species list and defvar_print_list
  # species_list: only species
  # defvar_print_list: what should be printed, maybe including some lines like
  #                    "{ Peroxy radicals. }"
  species_list = []
  defvar_print_list = []
  
  # Get the defvar species from the root def file
  update_defvar_list(file_lines, species_list, defvar_print_list)
  
  # Get the defvar species from the included files
  for file_name in include_file_list:
    with open(file_name, 'r') as f:
      tmp_file_lines = f.readlines()
      update_defvar_list(tmp_file_lines, species_list, defvar_print_list)
  
  #
  # Replace the #DEFVAR part in the root def file with new defvar_print_list
  #
  line_number_defvar, i1, i2 = index_containing_substring(file_lines, r'^\s*#DEFVAR', False)
  
  # Delete the old #DEFVAR block
  block_line_count, i1, i2 = index_containing_substring( \
    file_lines[line_number_defvar+1:], r'^\s*#', False)
  if block_line_count >= 0:
    del file_lines[line_number_defvar+1:line_number_defvar+block_line_count+1]
  else: # '#DEFVAR' is the last block
    del file_lines[line_number_defvar+1:]
  
  # Insert the new #DEFVAR block
  file_lines[line_number_defvar+1:line_number_defvar+1] = defvar_print_list
  
  # Write to log file
  logging.info('')
  logging.info('Combined all the species defined in #DEFVAR.')
  
  
  #---------- 3. Combine chemical species defined for #DEFFIX ----------#
  
  # Initiate species list and deffix_print_list
  # species_list: only species
  # defvar_print_list: what should be printed, maybe including some lines like
  #                    "{ Peroxy radicals. }"
  species_list = []
  deffix_print_list = []
  
  # Get the defvar species from the root def file
  update_deffix_list(file_lines, species_list, deffix_print_list)
  
  # Get the defvar species from the included files
  for file_name in include_file_list:
    with open(file_name, 'r') as f:
      tmp_file_lines = f.readlines()
      update_deffix_list(tmp_file_lines, species_list, deffix_print_list)
  
  #
  # Replace the #DEFFIX part in the root def file with new deffix_print_list
  #
  line_number_deffix, i1, i2 = index_containing_substring(file_lines, r'^\s*#DEFFIX', False)
  
  # If #DEFFIX exists in the original def file
  if line_number_deffix >= 0:
    # Delete the old #DEFFIX block
    block_line_count, i1, i2 = index_containing_substring( \
      file_lines[line_number_deffix+1:], r'^\s*#', False)
    if block_line_count >= 0:
      del file_lines[line_number_deffix+1:line_number_deffix+block_line_count+1]
    else: # '#DEFFIX' is the last block
      del file_lines[line_number_deffix+1:]
    
    # Insert the new #DEFFIX block
    file_lines[line_number_deffix+1:line_number_deffix+1] = deffix_print_list
  # Else create a new #DEFFIX block before #DEFVAR
  else:
    # Find the #DEFVAR
    line_number_defvar, i1, i2 = index_containing_substring(file_lines, r'^\s*#DEFVAR', False)

    # Insert the new #DEFFIX block before #DEFVAR
    deffix_print_list.insert(0, '#DEFFIX \n')
    file_lines[line_number_defvar-1:line_number_defvar-1] = deffix_print_list
  
  # Write to log file
  logging.info('')
  logging.info('Combined all the species defined in #DEFFIX.')
  
  
  #---------- 4. Combine RO2 species ----------#
  
  # Initiate RO2 list
  ro2_list = []
  
  # Update RO2 species for the root def file
  update_ro2_list(file_lines, ro2_list)
  
  # Get the RO2 species from the included files
  for file_name in include_file_list:
    with open(file_name, 'r') as f:
      tmp_file_lines = f.readlines()
      update_ro2_list(tmp_file_lines, ro2_list)
  
  # Write the new RO2 definition block according to new ro2_list
  ncol = 4
  nro2 = len(ro2_list)
  ro2_print_string = ''
  for i, ro2 in enumerate(ro2_list):
    # First one in a line
    if i % ncol == 0:
      ro2_print_string += '    C(ind_{0})'.format(ro2)
    # Others
    else:
      ro2_print_string += 'C(ind_{0})'.format(ro2)
  
    # Not the last RO2
    if i < nro2-1:
      ro2_print_string += ' + '
    
      # Last one in a line
      if (i+1) % ncol == 0:
        ro2_print_string += '&\n'
    # Last RO2
    else:
      ro2_print_string += '\n'
  
  #
  # Replace the RO2 definition in the root def file with new ro2_list
  #
  line_number_ro2, i1, i2 = index_containing_substring(file_lines, r'^\s*RO2\s*=', False)
  
  # Count the old RO2 definition block lines
  block_line_count = 0
  for i, s in enumerate(file_lines[line_number_ro2+1:]):
    # Count empty lines or the lines with 'C(ind_xxx)'
    if re.search(r'^\s*$', s) or re.search(r'C\(ind_', s, re.I):
      block_line_count += 1
    # Otherwise quit
    else:
      break
  
  # Delete the old RO2 definition block
  if block_line_count >= 0:
    del file_lines[line_number_ro2+1:line_number_ro2+block_line_count+1]
  
  # Write the new RO2 list
  file_lines[line_number_ro2+1:line_number_ro2+1] = ro2_print_string
  
  # Write to log file
  logging.info('')
  logging.info('Combined all the RO2 species.')
  
  
  #---------- 5. Combine chemical reactions for #EQUATIONS ----------#
  reactant_list = []
  product_list = []
  ratecoef_list = []
  equation_list = []
  equation_norate_list = []
  equation_print_list = []
  equation_duplicate_print_list = []
  
  # Update equations lists for the root def file
  update_equation_list(file_lines, equation_list, equation_norate_list, equation_print_list, equation_duplicate_print_list)
  
  # Update equations lists for the included files
  for file_name in include_file_list:
    with open(file_name, 'r') as f:
      tmp_file_lines = f.readlines()
      update_equation_list(tmp_file_lines, equation_list, equation_norate_list, equation_print_list, equation_duplicate_print_list)
  
  # Count the old #EQUATIONS block lines
  line_number_equations, i1, i2 = index_containing_substring(file_lines, r'^\s*#EQUATIONS', False)
  
  block_line_count_equations = 0
  for i, s in enumerate(file_lines[line_number_equations+1:]):
    if not re.search(r'^\s*#', s):  # not meet other blocks
      block_line_count_equations += 1
  
  # Delete the old #EQUATIONS block
  if block_line_count_equations >= 0:
    del file_lines[line_number_equations+1:line_number_equations+block_line_count_equations+1]
  
  # Write the new #EQUATIONS block
  file_lines[line_number_equations+1:line_number_equations+1] = equation_print_list
  
  # Write the duplicated equations to the log file
  logging.info('')
  logging.info('Combined all the equations:')
  logging.info('If the reactants, products and k values of several equations are the same, only one is left.')
  logging.info('If the reactants, products of several equations are the same, but k values are different, all of them are kept,')
  logging.info('and the followed duplicated equations will be shown below.')
  logging.info('The duplicated equations are:')
  for e in equation_duplicate_print_list:
    logging.info(re.sub(r'\n$', '', e))
  
  
  #==============================================================================#
  #
  # Update reaction rates of Criegee intermediate (CI) radicals producing SO3
  #
  # For example:
  # APINBOO + SO2 = PINAL + SO3 :  7.00D-14  ;
  # 7.00D-14 --> 7.00D-13
  #
  #==============================================================================#
  if is_ciso3:
    # Write to log file
    logging.info('')
    logging.info('The reaction rates of these equations (CI -> SO3) are modified from 7.00D-14 to 7.00D-13:')
  
    for i, s in enumerate(file_lines):
      if re.search(r'.*=.*SO3\s*:\s*7.00D-14', s, re.I):
        # Replace
        file_lines[i] = re.sub('7.00D-14', '7.00D-13', s)
  
        # Write to the log
        logging.info(s)
  
  
  #==============================================================================#
  #
  # Write the new contents to a new file
  #
  #==============================================================================#

  #---------- Remove the return characters ^M ----------#

  # The line ending '\r', '\n', '\r\n' will be translated automatically to '\n'
  # '\r' is just '^M'
  file_lines = [l.replace('\r', '') for l in file_lines]

  # Write to log file
  logging.info('')
  logging.info('Removed all the carriage return character ^M.')
  with open(args.modified_mcm_kpp_file_name, 'w') as f:
    f.writelines(file_lines)
  
  # Write to log file
  logging.info('')
  logging.info('The final file is written.')
