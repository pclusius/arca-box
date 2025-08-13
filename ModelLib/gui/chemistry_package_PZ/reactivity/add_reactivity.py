#!/usr/bin/python
"""
=============================================================================
Copyright (C) 2021  Multi-Scale Modelling group
Written by Zhou Putian, University of Helsinki
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

#==============================================================================#
#
#  Generate a second_reactivity.f90 file as a default name from the KPP
#  chemistry file.
#
#  Contact: Zhou Putian, putian.zhou@helsinki.fi, 19.11.2020
#
#  In the reactivity configuration file, every reactivity name is specified with
#  a species (e.g., OH, NO3, O3, etc.) and related reactants. There are also 2
#  special names for chemicals: +other and +all. +all will contain all chemicals
#  that react with the species in the chemistry file. +other will contain the
#  chemicals from +all that are not mentioned under any other reactivities for
#  this species. +all and +other must be the only chamical name in their class.
#
#==============================================================================#


#==============================================================================#
#
# Header
#
#==============================================================================#
import sys
import re  # regural expressions

import logging
import argparse
from argparse import RawTextHelpFormatter  # can insert new line in the text

import xml.etree.ElementTree as ET


#==============================================================================#
#
# Functions
#
#==============================================================================#
def find_reactions(kppdef_file, species):
  """
  " Return a dict containing the reactant names and corresponding reaction rate
  " list.
  " name: [rate1, rate2, ...]
  """

  species_dict = dict([])

  file = open(kppdef_file, 'r')

  # Forward the file until #EQUATIONS line
  for line in file:
    if re.search('^#EQUATIONS',line): break

  # Search equation lines
  for line in file:
    if re.search(r'^\s*$',line): continue  # empty line
    if re.search(r'^//',line): continue    # comment line

    equation_str = r'(^|\s|\+){0}(\s|\+|(?==)).*='.format(species)
    if re.search(equation_str,line):
      if re.search(r'\+.*\+.*=',line):
        logging.error('This program can only handle equations with species' \
            + ' and one other term on\n' \
            + 'the left side. Not more than one, like in this line:')
        logging.error(line)
        sys.exit()

    # Find the reactants with species
    if re.search('\{', line):  # with {tag} in the beginning
      equation_str = r'^\s*\{{.*?\}}\s*(?:{0}\s*\+\s*(\S+)|(\S+)\s*\+\s*{0})\s*=.*?:\s*?(.*)\s*?;'.format(species)
      m = re.search(equation_str, line)
    else:
      equation_str = r'^\s*(?:{0}\s*\+\s*(\S+)|(\S+)\s*\+\s*{0})\s*=.*?:\s*?(.*)\s*?;'.format(species)
      m = re.search(equation_str, line)

    # thanks to diogenes (Heikki Hallamaa) and wolverian (Ilmari Vacklin)
    # for help with this
    if m:
      chemical = m.groups()[0] or m.groups()[1]
      rate = m.groups()[2].strip()
      if chemical in species_dict:
        species_dict[chemical].append(rate)
      else:
        species_dict[chemical] = [rate]
  return species_dict


def read_reactivity_config(config_file):
  # Initiate lists
  tags = []
  species = []
  reactants = []

  # Get the root of xml file
  root = ET.parse(config_file).getroot()

  for reactivity in root.iter('reactivity'):
    # Get the reactivity names as tags
    tags.append(reactivity.find('name').text.strip())

    # Get the species name, e.g., OH, NO3, O3
    species.append(reactivity.find('species').text.strip())

    # Get the reactants list for each tag
    reactants_string = reactivity.find('reactants').text
    reactants_list = [i.strip() for i in reactants_string.split(',')]
    reactants.append(reactants_list)

  return tags, species, reactants


def create_second_reactivity(file_name, tags, reactants):
  # Variable definition list
  # var_def_list = ', '.join(tags)
  ntags = len(tags)

  with open(file_name, 'w') as f:
    f.write('module second_reactivity\n')
    f.write('\n')
    f.write('  implicit none\n')
    f.write('\n')
    f.write('  private\n')
    f.write('\n')
    f.write('  integer, parameter :: NREACTIVITY={0}\n'.format(ntags))
    f.write('\n')

    # Save all the reactivity names
    f.write('  character(len=20), parameter :: reactivity_name(NREACTIVITY) = (/ &\n')
    for i, t in enumerate(tags):
      if i==len(tags)-1:
        f.write("    '{0:20s}' &\n".format(t))
      else:
        f.write("    '{0:20s}', &\n".format(t))

      # # One name per line
      # if i == len(tags)-1:
      #   f.write(' /)\n')
      # else:
      #   f.write(', &\n')
    if len(tags) == 0: f.write('character(len=20):: &\n')
    f.write(' /)\n')

    # Public ones
    f.write('\n')
    f.write('  public :: calculate_reactivities, NREACTIVITY, reactivity_name\n')
    f.write('\n')
    f.write('contains\n')
    f.write('\n')

    # Begin the subroutine to calculate the reactivities
    f.write('  subroutine calculate_reactivities(CONC, reactivity)\n')

    # Use variables declared in other modules
    f.write('\n')
    f.write('    use second_Precision   ! dp\n')
    f.write('    use second_Parameters  ! NSPEC, ind_*\n')
    f.write('    use second_Global      ! TEMP, M, O2, N2, H2O, RO2, K values\n')

    # Declaration
    f.write('\n')
    f.write('    real(dp), intent(in) :: CONC(NSPEC)\n')
    f.write('    real(dp), intent(out) :: reactivity(NREACTIVITY)\n')
    f.write('\n')

    # Add the reaction rates together
    ncol = 5  # print 4 reactants in one line
    for it, (tag, s, reactant_list) in enumerate(zip(tags, species, reactants)):
      reaction_rate_lines = []
      reaction_rate_lines.append('    reactivity({0}) ='.format(it+1))

      for i, reactant in enumerate(reactant_list):
        if not reactant in species_dict[s]:
          logging.error('Chemical {0} (class: {1}) not found in any reaction equation'.format(reactant, tag))
          sys.exit()

        # Add all the reaction rates for this reactant
        total_rate = ' + '.join(species_dict[s][reactant])

        # Multiply the total rate with the reactant concentration
        reactivity_str = 'CONC(ind_{0})*({1})'.format(reactant, total_rate)

        # Do not put + before the first reactant
        if i == 0:
          reactivity_str = ' ' + reactivity_str
        else:
          reactivity_str = ' + ' + reactivity_str

        # Change a new line every 5 reactants
        if i % ncol == 0:
          reactivity_str = ' &\n     ' + reactivity_str

        # Add this reactant reactivity calculation string to the lines
        reaction_rate_lines.append(reactivity_str)

      # Write comments for this tag
      f.write('    ! {0}\n'.format(tag))

      # Write all the reactivities for this tag
      f.writelines(reaction_rate_lines)
      f.write('\n\n')

    # End the subroutine to calculate the reactivities
    f.write('  end subroutine calculate_reactivities\n')
    f.write('\n')

    # End of the module
    f.write('end module second_reactivity')


######## MAIN PROGRAM ########
if __name__ == '__main__':


  #===== Define the argument parser =====#
  help_string = \
    'Creating a second reactivity file to calculate the reactivities.\n' + \
    '\n' + \
    '1. Define the reactivity name, species and specify the related reactants.\n' + \
    'in the <config_file> xml file.\n' + \
    '\n' + \
    '2. Run this script with <config_file> to produce a Fortran code\n' + \
    '<second_file> which includes the declaration of number of reactivities\n' + \
    '(NREACTIVITY), reactivity name array (reactivity_name), and a subroutine\n' + \
    'calculate_reactivities(CONC, reactivity).\n' + \
    '\n' + \
    '3. Put the <second_file> to the chemistry scheme folder together with\n' + \
    'other second files, modify your Makefile when necessary.\n' + \
    '\n' + \
    '4. You can call the subroutine calculate_reactivities anywhere you want,\n' + \
    'but the suggested calling place is after the subroutine CHEMISTRY or\n' + \
    'after KPP_Proceed inside CHEMISTRY. Because in this case the variables\n' + \
    'TEMP, CONC, reaction rates and so on are updated.'

  parser = argparse.ArgumentParser(description=help_string, \
      add_help=True, formatter_class=RawTextHelpFormatter)

  default_values = {}
  default_values['config_file'] = 'reactivity.xml'
  default_values['kppdef_file'] = 'second.def'
  default_values['second_file'] = 'second_reactivity.f90'
  default_values['log_file'] = 'second_reactivity.log'

  parser.add_argument('config_file', \
      action='store', \
      help='Reactivity config file in xml format.')
  parser.add_argument('-k', \
      dest='kppdef_file', \
      default=default_values['kppdef_file'], \
      action='store', \
      help='The kpp def file which contains all the equations.\n' + \
        'The default is "{0}".'.format(default_values['kppdef_file']))
  parser.add_argument('-o', \
      dest='second_file', \
      default=default_values['second_file'], \
      action='store', \
      help='Output second_reactivity file. The default is\n' + \
        '"{0}".'.format(default_values['second_file']))
  parser.add_argument('--log', '-l', \
      dest='log_file', \
      default=default_values['log_file'], \
      action='store', \
      help='Log file.' + \
        ' The default is "{0}".'.format(default_values['log_file']))

  args = parser.parse_args()

  #===== Set up logging to file =====#
  logging.basicConfig(level=logging.DEBUG, \
                      format='%(asctime)s-12s %(levelname)-8s %(message)s', \
                      datefmt='%m-%d %H:%M',
                      filename=args.log_file,
                      filemode='w')

  # define a Handler which writes INFO messages or higher to the sys.stderr
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)

  # set a format which is simpler for console use
  formatter = logging.Formatter('%(asctime)-12s: %(levelname)-8s %(message)s')

  # tell the handler to use this format
  console.setFormatter(formatter)

  # add the handler to the root logger
  logging.getLogger().addHandler(console)

  # Output the input to the screen
  logging.info('  reactivity config file       : {0}'.format(args.config_file))
  logging.info('  input kpp def file           : {0}'.format(args.kppdef_file))
  logging.info('  output second_reactivity file: {0}'.format(args.second_file))

  #===== Read the config file =====#

  # Save reactivity names to tags, and other information
  tags, species, reactants = read_reactivity_config(args.config_file)

  # Get all the reactants and their reaction rates for a species
  species_dict = {}
  for s in list(set(species)):
    species_dict[s] = find_reactions(args.kppdef_file, s)

  # Get all the reactants for a species
  all_reactants_dict = {}
  for s in list(set(species)):
    all_reactants_dict[s] = list(species_dict[s].keys())

  # Find all the reactants for a species mentioned in the config file
  # sum(list, []): a workaround of flattening a list
  defined_reactants_dict = {}

  for s in list(set(species)):
    defined_reactants_dict[s] = []

  for s, r in zip(species, reactants):
    defined_reactants_dict[s].extend(sum(reactants, []))

  # Find all the other reactants for a species not mentioned in the config file
  other_reactants_dict = {}
  for s in list(set(species)):
    all_unique_set = set(all_reactants_dict[s])
    defined_unique_set = set(defined_reactants_dict[s])
    other_reactants_dict[s] = list(all_unique_set.difference(defined_unique_set))

  # Set the reactants list for '+all' and '+other'
  tmp_reactants = list(reactants)
  for i, (s, r) in enumerate(zip(species, reactants)):
    if r == ['+other']:
      tmp_reactants[i] = other_reactants_dict[s]
    if r == ['+all']:
      tmp_reactants[i] = all_reactants_dict[s]
  reactants = list(tmp_reactants)

  for i in range(len(tags)):
    logging.info('  class: {0} found {1} chemicals'.format(tags[i], len(reactants[i])))

  #===== Write out the second_file =====#
  create_second_reactivity(args.second_file, tags, reactants)
