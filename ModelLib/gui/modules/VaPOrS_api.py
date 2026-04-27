
"""
SMILES-Based Functional Group Identification and Vapor Pressure Calculation
=========================================================================

Created on: April 12, 2024
Author: Mojtaba Bezaatpour

### Overview
This script processes SMILES notations of organic compounds from an external file specified by the `filename_in` variable.
It identifies and quantifies 30 structural groups based on the SIMPOL method (DOI: 10.5194/acp-8-2773-2008),
calculates temperature-dependent saturation pressure and enthalpy of vaporization, and exports the results to an Excel file.

### Functional Groups Identified
The script detects and counts the following 30 structural groups:
- Molecular carbon
- Alkyl hydroxyl, aromatic hydroxyl
- Alkyl ether, alkyl ring ether, aromatic ether
- Aldehyde, ketone, carboxylic acid, ester
- Nitrate, nitro
- Alkyl amines (primary, secondary, tertiary), aromatic amine
- Amides (primary, secondary, tertiary)
- Peroxide, hydroperoxide, peroxy acid
- Non-aromatic C=C, carbonylperoxy acid
- Nitrate, nitro-phenol, nitro-ester
- Aromatic rings, non-aromatic rings
- C=C–C=O in a non-aromatic ring
- Carbon on the acid-side of an amide

### Functionality
- Extracts SMILES notations from an input file (`filename_in`).
- Identifies and quantifies the functional groups using the SIMPOL method.
- Computes temperature-dependent saturation pressure and enthalpy of vaporization for the compounds.
- Calculates vapor pressure and enthalpy of vaporization at a user-defined temperature (T).

### Output
The processed data is saved as `output.csv`, containing the identified functional groups and calculated vapor pressures.

"""



import numpy as np
from scipy.optimize import curve_fit

def find_closing_parenthesis(s, open_index):
    count = 0

    for i in range(open_index, len(s)):
        if s[i] == '(':
            count += 1
        elif s[i] == ')':
            count -= 1

            if count == 0:
                return i  # Found the corresponding closing parenthesis

    return -1  # Corresponding closing parenthesis not found


def find_opening_parenthesis(s, close_index):
    count = 0

    for i in range(close_index, -1, -1):
        if s[i] == ')':
            count += 1
        elif s[i] == '(':
            count -= 1

        if count == 0:
            return i  # Found the corresponding opening parenthesis

    return -1  # Corresponding opening parenthesis not found


def find_highest_digit(s):
    current_digit = 0
    highest_digit = None

    for char in s:
        if char.isdigit():
            current_digit = current_digit * 10 + int(char)
        elif current_digit != 0:
            if highest_digit is None or current_digit > highest_digit:
                highest_digit = current_digit
            current_digit = 0

    # Check for any remaining number at the end of the string
    if current_digit != 0:
        if highest_digit is None or current_digit > highest_digit:
            highest_digit = current_digit

    return highest_digit


def find_cycle_number(s):
    current_number = 0
    cycle_number = None
    number_count = {}  # Dictionary to keep track of number occurrences
    largest_single_digit = None  # Variable to track the largest single-digit number

    for char in s:
        if char.isdigit():
            current_number = current_number * 10 + int(char)
        elif current_number != 0:
            # Update the count of the current number
            if current_number in number_count:
                number_count[current_number] += 1
            else:
                number_count[current_number] = 1

            # Reset current_number for the next potential number
            current_number = 0

        # Check if the character is a single digit
        if char.isdigit() and len(char) == 1:
            digit = int(char)
            if largest_single_digit is None or digit > largest_single_digit:
                largest_single_digit = digit

    # Check for any remaining number at the end of the string
    if current_number != 0:
        if current_number in number_count:
            number_count[current_number] += 1
        else:
            number_count[current_number] = 1

    # Determine the cycle number
    for number, count in number_count.items():
        if count >= 2:
            # If found at least twice, consider it for cycle_number
            if cycle_number is None or number > cycle_number:
                cycle_number = number

    # If no number appeared twice, use the largest single-digit number
    if cycle_number is None:
        cycle_number = largest_single_digit

    return cycle_number

###################   Carbon number  #######################################################
def carbon_number(s):
    count_C = s.count('C')
    count_c = s.count('c')

    carbon = count_C + count_c

    return carbon

###################   Carbon number on acid side of amide  #######################################################
def ASA_carbon_number(s):
    primary_amide_number = 0
    # primary amide as first characters of the SMILES:
    if s[0:6] == 'O=C(N)':
        primary_amide_number = primary_amide_number + 1

    if s[0:6] == 'NC(=O)':
        primary_amide_number = primary_amide_number + 1

    # primary amide as last characters of the SMILES:
    if s[-6:] == 'C(=O)N':
        primary_amide_number = primary_amide_number + 1

    elif s[-6:] == 'C(N)=O':
        primary_amide_number = primary_amide_number + 1

    if s[-1] == 'N':
        if s[-4:-1] == 'O=C':
            primary_amide_number = primary_amide_number + 1
        elif s[-2] == ')':
            index = s.rfind(')', 0, -1)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-3:opening_parenthesis] == 'O=C':
                primary_amide_number = primary_amide_number + 1

    if s[-2:] == '=O':
        if s[-4:-2] == 'NC':
            primary_amide_number = primary_amide_number + 1
        elif s[-3] == ')':
            index = s.rfind(')', 0, -2)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-2:opening_parenthesis] == 'NC':
                primary_amide_number = primary_amide_number + 1

    # primary amide as last characters of a branch in the SMILES:
    co_index = s.find('C(=O)N)', 0)
    while co_index != -1:
        primary_amide_number = primary_amide_number + 1
        co_index = s.find('C(=O)N)', co_index + 1)

    co_index = s.find('C(N)=O)', 0)
    while co_index != -1:
        primary_amide_number = primary_amide_number + 1
        co_index = s.find('C(N)=O)', co_index + 1)

    co_index = s.find('N)', 0)
    while co_index != -1:
        if s[co_index-3:co_index] == 'O=C':
            primary_amide_number = primary_amide_number + 1
        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-3:opening_parenthesis] == 'O=C':
                primary_amide_number = primary_amide_number + 1

        co_index = s.find('N)', co_index + 1)

    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-2:co_index] == 'NC':
            primary_amide_number = primary_amide_number + 1
        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-2:opening_parenthesis] == 'NC':
                primary_amide_number = primary_amide_number + 1
        co_index = s.find('=O)', co_index + 1)
    if primary_amide_number == 1:
        count_C = s.count('C')
        count_c = s.count('c')
        asa_carbon_number = count_C + count_c
    else:
        asa_carbon_number = 0


    secondary_amide_number = 0
    # secondary amide as first characters of the SMILES:
    if s[0:3] == 'O=C':
        if s[3] == 'N':
            if s[4:6] == 'C(':
                secondary_amide_number = secondary_amide_number + 1
                asa_carbon_number = 1
        elif s[3] == '(':
            if s[4:6] == 'NC':
                secondary_amide_number = secondary_amide_number + 1
                index = s.find('(', 3)
                last_closing_parenthesis = find_closing_parenthesis(s, index)
                substring = s[last_closing_parenthesis:]
                asa_carbon_number = substring.count('C') + substring.count('c') + 1

            else:
                index = s.find('(', 3)
                closing_parenthesis = find_closing_parenthesis(s, index)
                if s[closing_parenthesis+1:closing_parenthesis+3] == 'NC':
                    secondary_amide_number = secondary_amide_number + 1
                    substring = s[3:closing_parenthesis]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

    # secondary amide as last characters of the SMILES:
    if s[-2:] == '=O':
        if s[-4:-2] == 'NC':
            if s[-5] == 'C':
                secondary_amide_number = secondary_amide_number + 1
                asa_carbon_number = 1
            elif s[-5] == ')':
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':
                    secondary_amide_number = secondary_amide_number + 1
                    asa_carbon_number = 1
                elif s[first_opening_parenthesis-1] == ')':
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':
                        secondary_amide_number = secondary_amide_number + 1
                        asa_carbon_number = 1

        elif s[-3] == ')':
            index = s.rfind(')', 0, -2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1:first_opening_parenthesis+3] == 'C(NC':
                secondary_amide_number = secondary_amide_number + 1
                substring = s[:first_opening_parenthesis]
                asa_carbon_number = substring.count('C') + substring.count('c')

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] == 'NC':
                if s[first_opening_parenthesis-3] == 'C':
                    secondary_amide_number = secondary_amide_number + 1
                    substring = s[first_opening_parenthesis:]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

                elif s[first_opening_parenthesis-3] == ')':
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':
                        secondary_amide_number = secondary_amide_number + 1
                        substring = s[first_opening_parenthesis:]
                        asa_carbon_number = substring.count('C') + substring.count('c') + 1

                    elif s[second_opening_parenthesis-1] == ')':
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':
                            secondary_amide_number = secondary_amide_number + 1
                            substring = s[first_opening_parenthesis:]
                            asa_carbon_number = substring.count('C') + substring.count('c') + 1

    # secondary amide as middle characters of the SMILES:
    co_index = s.find('NC(=O)', 0)
    while co_index != -1:
        if s[co_index-1] == 'C':
            secondary_amide_number = secondary_amide_number + 1
            substring = s[co_index:]
            asa_carbon_number = substring.count('C') + substring.count('c')

        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':
                secondary_amide_number = secondary_amide_number + 1
                substring = s[co_index:]
                asa_carbon_number = substring.count('C') + substring.count('c')

            elif s[first_opening_parenthesis-1] == ')':
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':
                    secondary_amide_number = secondary_amide_number + 1
                    substring = s[co_index:]
                    asa_carbon_number = substring.count('C') + substring.count('c')

        co_index = s.find('NC(=O)', co_index + 1)

    co_index = s.find('C(=O)NC', 0)
    while co_index != -1:
        secondary_amide_number = secondary_amide_number + 1
        if s[co_index-1] == '(':
            index = s.find('(', co_index-1)
            last_closing_parenthesis = find_closing_parenthesis(s, index)
            substring = s[:co_index] + s[last_closing_parenthesis:]
        else:
            substring = s[:co_index]

        asa_carbon_number = substring.count('C') + substring.count('c') + 1

        co_index = s.find('C(=O)NC', co_index + 1)

    co_index = s.find('(NC', 0)
    while co_index != -1:
        if s[co_index-1] == 'C':
            if s[co_index+3:co_index+6] == '=O)':
                secondary_amide_number = secondary_amide_number + 1
                asa_carbon_number = 1

            elif s[co_index+3:co_index+7] == '(=O)':
                secondary_amide_number = secondary_amide_number + 1
                index = s.find('(', co_index)
                last_closing_parenthesis = find_closing_parenthesis(s, index)
                substring = s[co_index:last_closing_parenthesis]
                asa_carbon_number = substring.count('C') + substring.count('c')

            elif s[co_index+3] == '(':
                index = s.find('(', co_index+3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':
                    secondary_amide_number = secondary_amide_number + 1
                    substring = s[co_index:first_closing_parenthesis]
                    asa_carbon_number = substring.count('C') + substring.count('c')

        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':
                if s[co_index+3:co_index+6] == '=O)':
                    secondary_amide_number = secondary_amide_number + 1
                    asa_carbon_number = 1

                elif s[co_index+3:co_index+7] == '(=O)':
                    secondary_amide_number = secondary_amide_number + 1
                    index = s.find('(', co_index)
                    last_closing_parenthesis = find_closing_parenthesis(s, index)
                    substring = s[co_index:last_closing_parenthesis]
                    asa_carbon_number = substring.count('C') + substring.count('c')

                elif s[co_index+3] == '(':
                    index = s.find('(', co_index+3)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':
                        secondary_amide_number = secondary_amide_number + 1
                        substring = s[co_index:first_closing_parenthesis]
                        asa_carbon_number = substring.count('C') + substring.count('c')

        co_index = s.find('(NC', co_index + 1)

    # secondary amide as last characters of a branch in the SMILES:
    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-2:co_index] == 'NC':
            if s[co_index-3] == 'C':
                secondary_amide_number = secondary_amide_number + 1
                asa_carbon_number = 1

            elif s[co_index-3] == ')':
                index = s.rfind(')', 0, co_index-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':
                    secondary_amide_number = secondary_amide_number + 1
                    asa_carbon_number = 1

                elif s[first_opening_parenthesis-1] == ')':
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':
                        secondary_amide_number = secondary_amide_number + 1
                        asa_carbon_number = 1

        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1:first_opening_parenthesis+3] == 'C(NC':
                secondary_amide_number = secondary_amide_number + 1
                substring = s[:first_opening_parenthesis]
                asa_carbon_number = substring.count('C') + substring.count('c')

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] == 'NC':
                if s[first_opening_parenthesis-3] == 'C':
                    secondary_amide_number = secondary_amide_number + 1
                    substring = s[first_opening_parenthesis:]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

                elif s[first_opening_parenthesis-3] == ')':
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':
                        secondary_amide_number = secondary_amide_number + 1
                        substring = s[first_opening_parenthesis:]
                        asa_carbon_number = substring.count('C') + substring.count('c') + 1

                    elif s[second_opening_parenthesis-1] == ')':
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':
                            secondary_amide_number = secondary_amide_number + 1
                            substring = s[first_opening_parenthesis:]
                            asa_carbon_number = substring.count('C') + substring.count('c') + 1

        co_index = s.find('=O)', co_index + 1)


    tertiary_amide_number = 0
    # tertiary amine as first characters of the SMILES:
    if s[0:7] == 'O=C(N(C':
        index = s.find('(', 5)
        closing_parenthesis = find_closing_parenthesis(s, index)
        if s[closing_parenthesis+1] == 'C':
            tertiary_amide_number = tertiary_amide_number + 1
            index0 = s.find('(', 3)
            last_closing_parenthesis = find_closing_parenthesis(s, index0)
            substring = s[last_closing_parenthesis:]
            asa_carbon_number = substring.count('C') + substring.count('c') + 1

    elif s[0:6] == 'O=CN(C':
        index = s.find('(', 4)
        closing_parenthesis = find_closing_parenthesis(s, index)
        if s[closing_parenthesis+1] == 'C':
            tertiary_amide_number = tertiary_amide_number + 1
            asa_carbon_number = 1

    elif s[0:4] == 'O=C(':
        index = s.find('(', 3)
        closing_parenthesis = find_closing_parenthesis(s, index)
        if s[closing_parenthesis+1:closing_parenthesis+4] == 'N(C':
            index2 = s.find('(', closing_parenthesis+2)
            second_closing_parenthesis = find_closing_parenthesis(s, index2)
            if s[second_closing_parenthesis+1] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                substring = s[index:closing_parenthesis]
                asa_carbon_number = substring.count('C') + substring.count('c') + 1

    # tertiary amine as last characters of the SMILES:
    if s[-4:] == ')C=O':
        index = s.rfind(')', 0, -3)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-1:first_opening_parenthesis+2] == 'N(C':
            if s[first_opening_parenthesis-2] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                asa_carbon_number = 1

            elif s[first_opening_parenthesis-2] == ')':
                index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                second_opening_parenthesis = find_opening_parenthesis(s, index1)
                if s[second_opening_parenthesis-1] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    asa_carbon_number = 1

                elif s[second_opening_parenthesis-1] == ')':
                    index2 = s.rfind(')', 0, second_opening_parenthesis)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        asa_carbon_number = 1

    if s[-3:] == ')=O':
        index = s.rfind(')', 0, -2)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-2:first_opening_parenthesis] == ')C':
            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
            second_opening_parenthesis = find_opening_parenthesis(s, index1)
            if s[second_opening_parenthesis-1:second_opening_parenthesis+2] == 'N(C':
                if s[second_opening_parenthesis-2] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    asa_carbon_number = 1

                elif s[second_opening_parenthesis-2] == ')':
                    index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        asa_carbon_number = 1

                    elif s[third_opening_parenthesis-1] == ')':
                        index3 = s.rfind(')', 0, third_opening_parenthesis)
                        forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                        if s[forth_opening_parenthesis-1] == 'C':
                            tertiary_amide_number = tertiary_amide_number + 1
                            asa_carbon_number = 1

        elif s[first_opening_parenthesis-1:first_opening_parenthesis+4] == 'C(N(C':
            index4 = s.find('(', first_opening_parenthesis+3)
            closing_parenthesis = find_closing_parenthesis(s, index4)
            if s[closing_parenthesis+1] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                substring = s[:first_opening_parenthesis-1]
                asa_carbon_number = substring.count('C') + substring.count('c') + 1

    # tertiary amine as middle characters of the SMILES:
    co_index = s.find(')C(=O)', 0)
    while co_index != -1:
        index = s.rfind(')', 0, co_index+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-1:second_opening_parenthesis+2] == 'N(C':
            if s[second_opening_parenthesis-2] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                substring = s[co_index+6:]
                asa_carbon_number = substring.count('C') + substring.count('c') + 1

            elif s[second_opening_parenthesis-2] == ')':
                index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                third_opening_parenthesis = find_opening_parenthesis(s, index2)
                if s[third_opening_parenthesis-1] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    substring = s[co_index+6:]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

                elif s[third_opening_parenthesis-1] == ')':
                    index3 = s.rfind(')', 0, third_opening_parenthesis)
                    forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                    if s[forth_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        substring = s[co_index+6:]
                        asa_carbon_number = substring.count('C') + substring.count('c') + 1

            elif s[second_opening_parenthesis-2] == '(':
                if s[second_opening_parenthesis-3] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    index = s.find('(', second_opening_parenthesis-2)
                    last_closing_parenthesis = find_closing_parenthesis(s, index)
                    substring = s[co_index+6:last_closing_parenthesis]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

                elif s[second_opening_parenthesis-3] == ')':
                    index2 = s.rfind(')', 0, second_opening_parenthesis-2)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        index = s.find('(', second_opening_parenthesis-2)
                        last_closing_parenthesis = find_closing_parenthesis(s, index)
                        substring = s[co_index+6:last_closing_parenthesis]
                        asa_carbon_number = substring.count('C') + substring.count('c') + 1

        co_index = s.find(')C(=O)', co_index + 1)

    co_index1 = s.find('C(=O)N(C', 0)
    while co_index1 != -1:
        index = s.find('(', 6)
        first_closing_parenthesis = find_closing_parenthesis(s, index)
        if s[first_closing_parenthesis+1] == 'C':
            tertiary_amide_number = tertiary_amide_number + 1
            if s[co_index1-1] == '(':
                index = s.find('(', co_index1-1)
                last_closing_parenthesis = find_closing_parenthesis(s, index)
                substring = s[:co_index1] + s[last_closing_parenthesis:]
            else:
                substring = s[:co_index1]

            asa_carbon_number = substring.count('C') + substring.count('c') + 1

        co_index1 = s.find('C(=O)N(C', co_index1 + 1)

    co_index2 = s.find(')C=O)', 0)
    while co_index2 != -1:
        index = s.rfind(')', 0, co_index2+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-2:second_opening_parenthesis+2] == '(N(C':
            if s[second_opening_parenthesis-3] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                asa_carbon_number = 1

            elif s[second_opening_parenthesis-3] == ')':
                index2 = s.rfind(')', 0, second_opening_parenthesis-2)
                third_opening_parenthesis = find_opening_parenthesis(s, index2)
                if s[third_opening_parenthesis-1] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    asa_carbon_number = 1

        co_index2 = s.find(')C=O)', co_index2 + 1)

    co_index3 = s.find(')=O)', 0)
    while co_index3 != -1:
        index = s.rfind(')', 0, co_index3+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-2:second_opening_parenthesis] == ')C':
            index2 = s.rfind(')', 0, second_opening_parenthesis-1)
            third_opening_parenthesis = find_opening_parenthesis(s, index2)
            if s[third_opening_parenthesis-2:third_opening_parenthesis+2] == '(N(C':
                if s[third_opening_parenthesis-3] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    substring = s[second_opening_parenthesis:co_index3]
                    asa_carbon_number = substring.count('C') + substring.count('c') + 1

                elif s[third_opening_parenthesis-3] == ')':
                    index3 = s.rfind(')', 0, third_opening_parenthesis-2)
                    forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                    if s[forth_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        substring = s[second_opening_parenthesis:co_index3]
                        asa_carbon_number = substring.count('C') + substring.count('c') + 1


        co_index3 = s.find(')=O)', co_index3 + 1)

    # tertiary amine as last characters of a branch in the SMILES:
    co_index4 = s.find(')C=O)', 0)
    while co_index4 != -1:
        index = s.rfind(')', 0, co_index4+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-1:first_opening_parenthesis+2] == 'N(C':
            if s[first_opening_parenthesis-2] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                asa_carbon_number = 1

            elif s[first_opening_parenthesis-2] == ')':
                index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                second_opening_parenthesis = find_opening_parenthesis(s, index1)
                if s[second_opening_parenthesis-1] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    asa_carbon_number = 1

                elif s[second_opening_parenthesis-1] == ')':
                    index2 = s.rfind(')', 0, second_opening_parenthesis)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        asa_carbon_number = 1

        co_index4 = s.find(')C=O)', co_index4 + 1)

    co_index5 = s.find(')=O)', 0)
    while co_index5 != -1:
        index = s.rfind(')', 0, co_index5+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-2:first_opening_parenthesis] == ')C':
            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
            second_opening_parenthesis = find_opening_parenthesis(s, index1)
            if s[second_opening_parenthesis-1:second_opening_parenthesis+2] == 'N(C':
                if s[second_opening_parenthesis-2] == 'C':
                    tertiary_amide_number = tertiary_amide_number + 1
                    asa_carbon_number = 1

                elif s[second_opening_parenthesis-2] == ')':
                    index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':
                        tertiary_amide_number = tertiary_amide_number + 1
                        asa_carbon_number = 1

                    elif s[third_opening_parenthesis-1] == ')':
                        index3 = s.rfind(')', 0, third_opening_parenthesis)
                        forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                        if s[forth_opening_parenthesis-1] == 'C':
                            tertiary_amide_number = tertiary_amide_number + 1
                            asa_carbon_number = 1

        elif s[first_opening_parenthesis-1:first_opening_parenthesis+4] == 'C(N(C':
            index4 = s.find('(', first_opening_parenthesis+3)
            closing_parenthesis = find_closing_parenthesis(s, index4)
            if s[closing_parenthesis+1] == 'C':
                tertiary_amide_number = tertiary_amide_number + 1
                substring = s[:first_opening_parenthesis-1] + s[co_index5+3:]
                asa_carbon_number = substring.count('C') + substring.count('c') + 1

        co_index5 = s.find(')=O)', co_index5 + 1)

    return asa_carbon_number

###################   Aromatic ring  #######################################################
def aromatic_ring(s):
    cyc_number = find_highest_digit(s)

    if cyc_number is not None:
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number + 1)]

        # Initialize the total sum
        total_sum = 0

        # Calculate the sum, counting each substring only once
        for substring in list_of_arom_rings:
            if substring in s:
                total_sum += 1  # Count only once if substring exists
        aroring = total_sum

    else:
        aroring = 0

    return aroring

###################   Nonaromatic ring  #######################################################
def non_aromatic_ring(s):

    # Find the number of cycles in the SMILES string
    cyc_number = find_cycle_number(s)

    # Define allowed non-aromatic atoms and numeric placeholders
    list_of_nonarom_atoms = ['C', 'N', 'O']
    digits = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

    if cyc_number is not None:
        # Create list of cycle indices based on cycle count
        list_of_index = [str(i) for i in range(1, cyc_number + 1)]

        # Initialize the total sum
        total_sum = 0

        # Loop through each cycle index to count non-aromatic rings
        for substring in list_of_index:
            if substring in s:
                A = s.find(substring, 0)    # First occurrence of the index
                B = s.find(substring, A+1)   # First occurrence of the index

                # Check for non-aromatic atoms around the cycle digits
                if s[A-1] in list_of_nonarom_atoms:    # C1, O1 or N1
                    if s[B-1] in list_of_nonarom_atoms:    # C1...O1
                        total_sum += 1

                    elif s[B-1] in digits:    # C2...O12
                        if s[B-2] in list_of_nonarom_atoms:
                            total_sum += 1

                elif s[A-1] in digits:  # C12...
                    if s[A-2] in list_of_nonarom_atoms:  # C12...
                        if s[B-1] in list_of_nonarom_atoms:  # C12...C2
                            total_sum += 1  #

                        elif s[B-1] in digits:  # C12...C32
                            if s[B-2] in list_of_nonarom_atoms:
                                total_sum += 1  #

        # Return the count of non-aromatic rings
        nonaroring = total_sum

    else:
        nonaroring = 0

    return nonaroring

###################   Nonaromatic double_bond carbons   #######################################################
def double_bound_nonaromatic_carbons(s):
    # Initialize the count of double-bonded non-aromatic carbons
    double_bound_nonaromatic_number = 0

    # Find the highest digit indicating cycle number in the SMILES string
    cyc_number = find_highest_digit(s)

    # Generate a list of non-aromatic ring markers if cycle count is found
    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_nonarom_rings = []

    # Find the index of each '=C' pattern in the SMILES string
    co_index = s.find('=C', 0)

    while co_index != -1:
        if s[co_index-1] == 'C':  # C=C
            double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

        elif co_index-2 >= 0 and s[co_index-2:co_index] in list_of_nonarom_rings:  # C1=C
            double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

        elif co_index-3 >= 0 and s[co_index-3:co_index] in list_of_nonarom_rings:  # C12=C
            double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

        elif s[co_index-1] == '(':  # (=C
            if s[co_index-2] == 'C':  # C(=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

            elif co_index-3 >= 0 and s[co_index-3:co_index-1] in list_of_nonarom_rings:  # C1(=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

            elif co_index-4 >= 0 and s[co_index-4:co_index-1] in list_of_nonarom_rings:  # C12(=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

        elif s[co_index-1] == ')':  # )=C
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  # C(...)=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_nonarom_rings:  # C1(...)=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_nonarom_rings:  # C12(...)=C
                double_bound_nonaromatic_number = double_bound_nonaromatic_number + 1

        # Move to the next occurrence of '=C'
        co_index = s.find('=C', co_index + 1)


    return double_bound_nonaromatic_number

###################   C=C-C=O  #######################################################
def nonaromatic_CCCO(s):

    CCCO = 0   # Counter for C=C-C=O occurrences


    cyc_number = find_highest_digit(s)  # Identify the highest cycle digit in SMILES notation

    if cyc_number != None: # if there is any cycle
        for i in range(1,cyc_number+1):
            A = s.find(str(i), 0)    # First occurrence of cycle number
            B = s.find(str(i), A+1)   #  Second occurrence of cycle number

            #  C=CC=O as first characters of the SMILES
            if s.startswith('O=C'+str(i)):
                if len(s) >= 6 and s[4:6] == 'C(':  #  /O=C1C(
                    first_closing_parenthesis = find_closing_parenthesis(s, 5)
                    if len(s) >= 8 and s[6:8] == '=C' and B < first_closing_parenthesis:   #  /O=C1C(=C...1)...
                        CCCO = CCCO + 1

                    if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=C' and first_closing_parenthesis < B:  #  /O=C1C(...)=C...1
                        CCCO = CCCO + 1

                if len(s) >= 7 and s[4:7] == 'C=C':  #  /O=C1C=C
                        CCCO = CCCO + 1


                if s[B-2:B] == '=C':  #  /O=C1...=C1
                    if s[B-3] == 'C':  #  /O=C1...C=C1
                        CCCO = CCCO + 1

                    if s[B-3] == '(':  #  /O=C1...(=C1
                        if s[B-4] == 'C':  #  /O=C1...C(=C1
                            CCCO = CCCO + 1

                    if s[B-3] == ')':  #  /O=C1...)=C1
                        index = s.rfind(')', 0, B-2)
                        first_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[first_opening_parenthesis-1] == 'C':  #  /O=C1...C(...)=C1
                            CCCO = CCCO + 1

            #  C=CC=O as middle characters of the SMILES
            co_index = s.find('C'+str(i)+'=C', 0)
            while co_index != -1:
                if len(s) > co_index+4 and s[co_index+4] == 'C':  # ...C1=CC
                    if len(s) > co_index+5 and s[co_index+5] == '(':  #  ...C1=CC(
                        if len(s) > co_index+7 and s[co_index+6:co_index+8] == '=O':  #  ...C1=CC(=O)...1...
                            if B > co_index+8:  #  ...C1=CC(=O)...1...
                                CCCO = CCCO + 1

                        first_closing_parenthesis = find_closing_parenthesis(s, co_index+5)
                        if B > co_index+5 and B < first_closing_parenthesis:  #  ...C1=CC(...1...)...
                            if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  #  ...C1=CC(...1...)=O
                                CCCO = CCCO + 1

                if len(s) > co_index+4 and s[co_index+4] == '(':  # ...C1=C(
                    first_closing_parenthesis = find_closing_parenthesis(s, co_index+4)
                    if len(s) > co_index+6 and s[co_index+5:co_index+7] == 'C(':  # ...C1=C(C(
                        second_closing_parenthesis = find_closing_parenthesis(s, co_index+6)
                        if len(s) > co_index+8 and s[co_index+7:co_index+9] == '=O':  # ...C1=C(C(=O)...
                            if B > second_closing_parenthesis and B < first_closing_parenthesis:  # ...C1=C(C(=O)... 1...)...
                                CCCO = CCCO + 1

                        if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  # ...C1=C(C(...)=O
                            if B > co_index+6 and B < second_closing_parenthesis:  # ...C1=C(C(...1...)=O
                                CCCO = CCCO + 1

                    if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == 'C(':  # ...C1=C(...)C(
                        second_closing_parenthesis = find_closing_parenthesis(s, first_closing_parenthesis+2)
                        if len(s) > first_closing_parenthesis+4 and s[first_closing_parenthesis+3:first_closing_parenthesis+5] == '=O':  # ...C1=C(...)C(=O)
                            if B > second_closing_parenthesis:  # ...C1=C(...)C(=O)...1...
                                CCCO = CCCO + 1

                        if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  # ...C1=C(...)C(...)=O
                            if B > first_closing_parenthesis+2 and B < second_closing_parenthesis:  # ...C1=C(...)C(...1...)=O
                                CCCO = CCCO + 1

                    if s[B-1:B+3] == 'C'+str(i)+'=O':  # ...C1=C(...)...C1=O  or  # ...C1=C(...C1=O )...
                        CCCO = CCCO + 1

                co_index = s.find('C'+str(i)+'=C', co_index + 1)


            co_index2 = s.find('C(=O)C', 0)
            while co_index2 != -1:
                if len(s) > co_index2+7 and s[co_index2+6:co_index2+8] == '=C':  # ...C(=O)C=C
                    if A < co_index2 and B > co_index2+7:  # 1...C(=O)C=C...1
                        CCCO = CCCO + 1

                if co_index2 + 6 < len(s) and s[co_index2+6] == '(':  # ...C(=O)C(
                    first_closing_parenthesis = find_closing_parenthesis(s, co_index2+6)
                    if len(s) > co_index2+8 and s[co_index2+7:co_index2+9] == '=C':  # ...C(=O)C(=C
                        if A < co_index2 and B > co_index2+8 and B < first_closing_parenthesis:  # 1...C(=O)C(=C...1...)...
                            CCCO = CCCO + 1

                    if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=C':  # ...C(=O)C(...)=C
                        if A < co_index2 and B > first_closing_parenthesis+2:  # 1...C(=O)C()=C...1......
                            CCCO = CCCO + 1

                co_index2 = s.find('C(=O)C', co_index2 + 1)


            co_index3 = s.find('C(', 0)
            while co_index3 != -1:
                first_closing_parenthesis = find_closing_parenthesis(s, co_index3+1)
                if s[co_index3+2] == 'C':  #  C(C
                    if len(s) > co_index3+4 and s[co_index3+3:co_index3+5] == '=C':  #  C(C=C
                        if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  #  C(C=C...)=O
                            if A < co_index3 and B > co_index3+4 and B < first_closing_parenthesis:  #  1...C(C=C...1...)=O
                                CCCO = CCCO + 1


                    if len(s) > co_index3+3 and s[co_index3+3] == '(':  #  C(C(
                        second_closing_parenthesis = find_closing_parenthesis(s, co_index3+3)
                        if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  #  C(C(...)...)=O
                            if len(s) > co_index3+5 and s[co_index3+4:co_index3+6] == '=C':  #  C(C(=C...)...)=O
                                if A < co_index3 and B > co_index3+5 and B < second_closing_parenthesis:  #  C(C(=C...1...)...)=O
                                    CCCO = CCCO + 1

                            if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=C':  #  C(C(...)=C...)=O
                                if A < co_index3 and B > second_closing_parenthesis+2 and B < first_closing_parenthesis:  #  C(C(...)=C...1...)=O
                                    CCCO = CCCO + 1

                if len(s) > co_index3+3 and s[co_index3+2:co_index3+4] == '=C':  #  C(=C
                    if len(s) > co_index3+4 and s[co_index3+4] == 'C':  #  C(=CC
                        if len(s) > co_index3+5 and s[co_index3+5] == '(':  #  C(=CC(
                            second_closing_parenthesis = find_closing_parenthesis(s, co_index3+5)
                            if len(s) > co_index3+7 and s[co_index3+6:co_index3+8] == '=O':  #  C(=CC(=O)...
                                if A < co_index3 and B > second_closing_parenthesis and B < first_closing_parenthesis:  # C(=CC(=O)...1...)
                                    CCCO = CCCO + 1

                            if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  #  C(=CC(...)=O
                                if A < co_index3 and B > co_index3+5 and B < second_closing_parenthesis:  #  C(=CC(...1...)=O
                                    CCCO = CCCO + 1

                        if len(s) > co_index3+5 and s[co_index3+5] == str(i):  #  C(=CC1
                            if len(s) > co_index3+7 and s[co_index3+6:co_index3+8] == '=O':  #  C(=CC1=O...)...
                                if A < co_index3:  #  1...C(=CC1=O...)...
                                    CCCO = CCCO + 1

                    if len(s) > co_index3+4 and s[co_index3+4] == '(':  #  C(=C(
                        second_closing_parenthesis = find_closing_parenthesis(s, co_index3+4)
                        if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == 'C(':  #  C(=C(...)C(
                            third_closing_parenthesis = find_closing_parenthesis(s, second_closing_parenthesis+2)
                            if len(s) > second_closing_parenthesis+4 and s[second_closing_parenthesis+3:second_closing_parenthesis+5] == '=O':  #  C(=C(...)C(=O)
                                if A < co_index3 and B > second_closing_parenthesis+5 and B < first_closing_parenthesis:  #  C(=C(...)C(=O)...1...)
                                    CCCO = CCCO + 1

                            if len(s) > third_closing_parenthesis+2 and s[third_closing_parenthesis+1:third_closing_parenthesis+3] == '=O':  #  C(=C(...)C(...)=O)...
                                if A < co_index3 and B > second_closing_parenthesis+2 and B < third_closing_parenthesis:
                                    CCCO = CCCO + 1

                        if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == 'C'+str(i):  #  C(=C(...)C1
                            if len(s) > second_closing_parenthesis+4 and s[second_closing_parenthesis+3:second_closing_parenthesis+5] == '=O':  #  C(=C(...)C1=O
                                if A < co_index3:  #  1..C(=C(...)C1=O
                                    CCCO = CCCO + 1

                        if len(s) > co_index3+5 and s[co_index3+5] == 'C':  #  C(=C(C
                            if len(s) > co_index3+6 and s[co_index3+6] == str(i):  #  C(=C(C1
                                if len(s) > co_index3+8 and s[co_index3+7:co_index3+9] == '=O':  #  C(=C(C1=O...)...)..
                                    if A < co_index3:  #  1...C(=C(C1=O...)...)...
                                        CCCO = CCCO + 1

                            if len(s) > co_index3+6 and s[co_index3+6] == '(':  #  C(=C(C(
                                third_closing_parenthesis = find_closing_parenthesis(s, co_index3+6)
                                if len(s) > co_index3+9 and s[co_index3+7:co_index3+10] == '=O':  #  C(=C(C(=O)
                                    if A < co_index3 and B > co_index3+9 and B < second_closing_parenthesis:
                                        CCCO = CCCO + 1

                                if len(s) > third_closing_parenthesis+2 and s[third_closing_parenthesis+1:third_closing_parenthesis+3] == '=O':  #  C(=C(C(...)=O)...)
                                    if A < co_index3 and B > co_index3+6 and B < third_closing_parenthesis:  #  C(=C(C(...1...)=O)...)
                                        CCCO = CCCO + 1


                if len(s) > first_closing_parenthesis+2 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=C':  #  C(...)=C
                    if len(s) > first_closing_parenthesis+3 and s[first_closing_parenthesis+3] == '(':  #  C(...)=C(
                        second_closing_parenthesis = find_closing_parenthesis(s, first_closing_parenthesis+3)
                        if len(s) > second_closing_parenthesis+1 and s[second_closing_parenthesis+1] == 'C':  #  C(...)=C(...)C
                            if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+2] == str(i):  #  C(...)=C(...)C1
                                if len(s) > second_closing_parenthesis+4 and s[second_closing_parenthesis+3:second_closing_parenthesis+5] == '=O':  #  C(...)=C(...)C1=O
                                    if A < co_index3:  #  1...C(...)=C(...)C1=O
                                        CCCO = CCCO + 1

                            if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+2] == '(':  #  C(...)=C(...)C(
                                third_closing_parenthesis = find_closing_parenthesis(s, second_closing_parenthesis+2)
                                if len(s) > second_closing_parenthesis+5 and s[second_closing_parenthesis+3:second_closing_parenthesis+6] == '=O)':  #  C(...)=C(...)C(=O)
                                    if A < co_index3 and B > second_closing_parenthesis+5:  #  1...C(...)=C(...)C(=O)...1
                                        CCCO = CCCO + 1

                                if len(s) > third_closing_parenthesis+2 and s[third_closing_parenthesis+1:third_closing_parenthesis+3] == '=O':  #  C(...)=C(...)C(...)=O
                                    if A < co_index3 and B > second_closing_parenthesis+2 and B < third_closing_parenthesis:  #  C(...)=C(...)C(...1...)=O
                                        CCCO = CCCO + 1

                        if len(s) > first_closing_parenthesis+4 and s[first_closing_parenthesis+4] == 'C':  #  C(...)=C(C
                            if len(s) > first_closing_parenthesis+5 and s[first_closing_parenthesis+5] == str(i):  #  C(...)=C(C1
                                if len(s) > first_closing_parenthesis+7 and s[first_closing_parenthesis+6:first_closing_parenthesis+8] == '=O':  #  C(...)=C(C1=O...)...
                                    if A < co_index3:  #  1...C(...)=C(C1=O...)...
                                        CCCO = CCCO + 1

                            if len(s) > first_closing_parenthesis+5 and s[first_closing_parenthesis+5] == '(':  #  C(...)=C(C(
                                third_closing_parenthesis = find_closing_parenthesis(s, first_closing_parenthesis+5)
                                if len(s) > first_closing_parenthesis+8 and s[first_closing_parenthesis+6:first_closing_parenthesis+9] == '=O)':  #  C(...)=C(C(=O)...
                                    if A < co_index3 and B > first_closing_parenthesis+8 and B < second_closing_parenthesis:  #  1...C(...)=C(C(=O)...1...)...
                                        CCCO = CCCO + 1

                                if len(s) > third_closing_parenthesis+2 and s[third_closing_parenthesis+1:third_closing_parenthesis+3] == '=O':  #  C(...)=C(C(...)=O
                                    if A < co_index3 and B > first_closing_parenthesis+5 and B < third_closing_parenthesis:  #  1...C(...)=C(C(...1...)=O
                                        CCCO = CCCO + 1




                    if len(s) > first_closing_parenthesis+3 and s[first_closing_parenthesis+3] == 'C':  #  C(...)=CC
                        if len(s) > first_closing_parenthesis+4 and s[first_closing_parenthesis+4] == str(i):  #  C(...)=CC1
                            if len(s) > first_closing_parenthesis+6 and s[first_closing_parenthesis+5:first_closing_parenthesis+7] == '=O':  #  C(...)=CC1=O
                                if A < co_index3:  #  1...C(...)=CC1=O
                                    CCCO = CCCO + 1

                        if len(s) > first_closing_parenthesis+4 and s[first_closing_parenthesis+4] == '(':  #  C(...)=CC(
                            second_closing_parenthesis = find_closing_parenthesis(s, first_closing_parenthesis+4)
                            if len(s) > first_closing_parenthesis+7 and s[first_closing_parenthesis+5:first_closing_parenthesis+8] == '=O)':  #  C(...)=CC(=O)
                                if A < co_index3 and B > first_closing_parenthesis+7:  #  1...C(...)=CC(=O)...1
                                    CCCO = CCCO + 1

                            if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  #  C(...)=CC(...)=O
                                if A < co_index3 and B > first_closing_parenthesis+4 and B < second_closing_parenthesis:  #  1...C(...)=CC(...1...)=O
                                    CCCO = CCCO + 1


                co_index3 = s.find('C(', co_index3 + 1)

            co_index4 = s.find('C=C', 0)
            while co_index4 != -1:
                if len(s) >= co_index4+4 and s[co_index4+3] == 'C':  #  C=CC
                    if len(s) >= co_index4+5 and s[co_index4+4] == str(i):  #  C=CC1
                        if len(s) >= co_index4+7 and s[co_index4+5:co_index4+7] == '=O':  #  C=CC1=O
                            if A < co_index4:  #  1...C=CC1=O
                                CCCO = CCCO + 1

                    if len(s) >= co_index4+5 and s[co_index4+4] == '(':  #  C=CC(
                        first_closing_parenthesis = find_closing_parenthesis(s, co_index4+4)
                        if len(s) >= first_closing_parenthesis+3 and s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  #  C=CC(...)=O
                            if A < co_index4 and B > co_index4+4 and B < first_closing_parenthesis:  #  C=CC(...1...)=O
                                CCCO = CCCO + 1

                        if len(s) >= co_index4+7 and s[co_index4+5:co_index4+7] == '=O':  #  C=CC(=O)...
                            if A < co_index4 and B > co_index4+7:  #  1...C=CC(=O)...1
                                CCCO = CCCO + 1

                if len(s) >= co_index4+4 and s[co_index4+3] == '(':  #  C=C(
                    first_closing_parenthesis = find_closing_parenthesis(s, co_index4+3)
                    if len(s) >= co_index4+5 and s[co_index4+4] == 'C':  #  C=C(C
                        if len(s) >= co_index4+6 and s[co_index4+5] == str(i):  #  C=C(C1
                            if len(s) >= co_index4+8 and s[co_index4+6:co_index4+8] == '=O':  #  C=C(C1
                                if A < co_index4:  #  1...C=C(C1
                                    CCCO = CCCO + 1

                        if len(s) >= co_index4+6 and s[co_index4+5] == '(':  #  C=C(C(
                            second_closing_parenthesis = find_closing_parenthesis(s, co_index4+5)
                            if len(s) >= second_closing_parenthesis+3 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  #  C=C(C(...)=O)
                                if A < co_index4 and B > co_index4+5 and B < second_closing_parenthesis:  #  C=C(C(...1...)=O)
                                    CCCO = CCCO + 1

                            if len(s) >= co_index4+8 and s[co_index4+6:co_index4+8] == '=O':  #  C=C(C(=O)
                                if A < co_index4 and B > co_index4+8 and B < first_closing_parenthesis:  #  1...C=C(C(=O)...1...)...
                                    CCCO = CCCO + 1

                    if len(s) >= first_closing_parenthesis+2 and s[first_closing_parenthesis+1] == 'C':  #  C=C(...)C
                        if len(s) >= first_closing_parenthesis+3 and s[first_closing_parenthesis+2] == str(i):  #  C=C(...)C1
                            if len(s) >= first_closing_parenthesis+5 and s[first_closing_parenthesis+3:first_closing_parenthesis+5] == '=O':  #  C=C(...)C1=O
                                if A < co_index4:  #  1...C=C(...)C1=O
                                    CCCO = CCCO + 1

                        if len(s) >= first_closing_parenthesis+3 and s[first_closing_parenthesis+2] == '(':  #  C=C(...)C(
                            second_closing_parenthesis = find_closing_parenthesis(s, first_closing_parenthesis+2)
                            if len(s) >= second_closing_parenthesis+5 and s[first_closing_parenthesis+3:first_closing_parenthesis+5] == '=O':  #  C=C(...)C(=O)...
                                if A < co_index4 and B > first_closing_parenthesis+5:  #  1...C=C(...)C(=O)...1...
                                    CCCO = CCCO + 1

                            if len(s) >= second_closing_parenthesis+3 and s[second_closing_parenthesis+1:second_closing_parenthesis+3] == '=O':  #  C=C(...)C(...)=O...
                                if A < co_index4 and B > first_closing_parenthesis+2 and B < second_closing_parenthesis:  #  1...C=C(...)C(...1...)=O...
                                    CCCO = CCCO + 1


                co_index4 = s.find('C=C', co_index4 + 1)


    if cyc_number != None: # if there is any cycle
        # Ending with cycle
        for i in range(1, cyc_number+1):
            end_str = 'C' + str(i) + '=O'
            if s.endswith(end_str):  # ...C1=O/
                index = s.find(str(i), 0)
                if s[index-1] == 'C' and s[index + len(str(i)): index + len(str(i)) + 2] == '=C':  # C1=C...C1=O/
                    CCCO = CCCO + 1

        for i in range(1, cyc_number+1):
            co_index = s.find('C' + str(i) + '=O)', 0) # ...C1=O)
            while co_index != -1:
                index = s.find(str(i), 0)
                if s[index-1] == 'C' and s[index + len(str(i)): index + len(str(i)) + 2] == '=C':  # C1=C...C1=O)
                    CCCO = CCCO + 1

                co_index = s.find('C' + str(i) + '=O)', co_index + 1)


    return CCCO

###################   Hydroxyl  #######################################################
def hydroxyl_group(s):
    hydroxyl_number = 0
    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_rings = []

    # hydroxyl as last characters of the SMILE:
    if s[-1] == 'O':  #  O/
        if s[-2] == 'C' or s[-2] == 'N':  #  CO/ or NO/
            if s[-4:-2] != 'O=':  #  /O=CO/ !
                hydroxyl_number = hydroxyl_number + 1

        elif s[-2] == ')':  #  )O/
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C' or s[first_opening_parenthesis-1] == 'N':  #  C(...)O/ or #  N(...)O/
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  C(=O)O/ !
                    if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(...)O/ !
                        hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  )(...)O/
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C' or s[second_opening_parenthesis-1] == 'N':  #  C(...)(...)O/ or #  N(...)(...)O/
                    hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  #  C1(...)O/
                hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  #  C12(...)O/
                hydroxyl_number = hydroxyl_number + 1

        elif s[-3:-1] in list_of_rings:  #  C1O/
            hydroxyl_number = hydroxyl_number + 1

        elif s[-4:-1] in list_of_rings:  #  C12O/
            hydroxyl_number = hydroxyl_number + 1


    #  hydroxyl as first characters of the SMILE:
    #  excluding if the string starts with hydroxyl i.e., /CO/ since they were included above
    if s.startswith('C('):  #  /C(
        if s[2:4] == 'O)':  #  /C(O)
            if s[4:6] == '=O':  #  /C(O)=O!
                pass

            elif s[4] == '(':  #  /C(O)(
                if s[5:7] == '=O':  #  /C(O)(=O)!
                    pass

                elif s[5:7] == 'O)':  #  /C(O)(O)
                    if s[7:9] != '=O':  #  /C(O)(O)=O !
                        hydroxyl_number = hydroxyl_number + 2

                else:
                    index = s.find('(', 4)
                    second_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O':  #  /C(O)(...)=O!
                        hydroxyl_number = hydroxyl_number + 1

            else:
                hydroxyl_number = hydroxyl_number + 1

        elif s[2:4] == '=O':  #  /C(=O) !
            pass

        else:
            index = s.find('(', 1)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '(O)':  #  /C(...)(O)...
                if s[first_closing_parenthesis+4:first_closing_parenthesis+6] != '=O':  #  /C(...)(O)=O/ !
                    hydroxyl_number = hydroxyl_number + 1

    if s.startswith('OC'):  #  /OC
        if len(s) > 3 and s[2:4] == '=O':  #  /OC=O!
            pass

        elif len(s) > 2 and s[2] == '(':  #  /OC(
            if s[3:5] == '=O':  #  /OC(=O)!
                pass

            else:
                index = s.find('(', 2)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O': #  /OC(...)=O!
                    pass

                else:
                    hydroxyl_number = hydroxyl_number + 1  # /OC(O)... can have two hydroxyl groups which is considered in next parts below.

        else:
            hydroxyl_number = hydroxyl_number + 1


           # nonaromatic ring-connected hydroxyl
    for j in list_of_rings:
        if s.startswith(j + '('):  #  /C1(
            if s[3:5] == 'O)':  #  /C1(O)
                if s[5] == '(':  #  /C1(O)(
                    if s[6:8] == 'O)':  #  /C1(O)(O)
                        hydroxyl_number = hydroxyl_number + 2

                    else:
                        hydroxyl_number = hydroxyl_number + 1

                else:
                    hydroxyl_number = hydroxyl_number + 1

            elif s[3:5] == '=O':  #  /C1(=O) !
                pass

            else:
                index = s.find('(', 2)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '(O)':  #  /C1(...)(O)...
                    hydroxyl_number = hydroxyl_number + 1

        if s.startswith('O' + j):  #  /OC1
                hydroxyl_number = hydroxyl_number  # This is coverd by /OC pattern too (no need for counting)


    #  hydroxyl as middle characters of the SMILE:
    co_index = s.find('(O)', 0)
    while co_index != -1:
        if s[co_index-1] == 'C':  #  C(O)
            if co_index-1 != 0:  #  /C(O)!
                if s[co_index-3:co_index-1] == 'O=':  #  O=C(O)!
                    pass

                elif s[co_index+3:co_index+5] == '=O':  #  ...C(O)=O!
                    pass

                else:
                    hydroxyl_number = hydroxyl_number + 1

        elif s[co_index-2:co_index] in list_of_rings:  #  C1(O)...
            if co_index-2 != 0:  #  /C1(O)...!
                hydroxyl_number = hydroxyl_number + 1

        elif s[co_index-3:co_index] in list_of_rings:  #  C12(O)...
            hydroxyl_number = hydroxyl_number + 1

        elif s[co_index-1] == ')':  #  )(O)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  C(...)(O)...
                if first_opening_parenthesis-1 != 0:  #  /C(...)(O)...!
                    hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  #  C1(...)(O)...
                if first_opening_parenthesis-2 != 0:  #  /C1(...)(O)...!
                    hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  #  C12(...)(O)...
                hydroxyl_number = hydroxyl_number + 1

        co_index = s.find('(O)', co_index + 1)


    #  hydroxyl as last characters of a branch in the SMILE:
    co_index2 = s.find('O)', 0)
    while co_index2 != -1:
        if s[co_index2-1] == 'C':  # CO)
            hydroxyl_number = hydroxyl_number + 1

        elif s[co_index2-2:co_index2] in list_of_rings:  # C1O)
            hydroxyl_number = hydroxyl_number + 1

        elif s[co_index2-3:co_index2] in list_of_rings:  # C12O)
            hydroxyl_number = hydroxyl_number + 1

        elif s[co_index2-1] == ')':  # )O)
            index = s.rfind(')', 0, co_index2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  # ...C(...)O)
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  # ...C(=O)O)!
                    hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  # ...C1(...)O)
                hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  # ...C12(...)O)
                hydroxyl_number = hydroxyl_number + 1

            elif s[first_opening_parenthesis-1] == ')':  # ...)(...)O)
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  # C(...)(...)O)
                    hydroxyl_number = hydroxyl_number + 1

        co_index2 = s.find('O)', co_index2 + 1)


    return hydroxyl_number

###################   Aldehyde  #######################################################
def aldehyde_group(s):
    aldehyde_number = 0

    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings

    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # aldehyde as first characters of the SMILE
    if s.startswith('O=CC') or s.startswith('C(=O)C') or s.startswith('O=Cc') or s.startswith('C(=O)c'):  #/O=CC or /C(=O)C or /O=Cc or /C(=O)c
        aldehyde_number = aldehyde_number + 1

    if s.startswith('O=C'):
        if len(s) == 3:
            aldehyde_number = aldehyde_number + 1

    # aldehyde as last characters of the SMILE
    if len(s) >= 3 and s[-3:] == 'C=O':  # C=O/
        if len(s) >= 4 and s[-4] == 'C':  # CC=O/
            aldehyde_number = aldehyde_number + 1

        # aldehyde connected to a ring
        elif s[-5:-3] in list_tot:  # C1C=O/ or c1C=O/
            aldehyde_number = aldehyde_number + 1

        elif s[-6:-3] in list_tot:  # C12C=O/ or c13C=O/
            aldehyde_number = aldehyde_number + 1

        elif s[-4:] == ')C=O':  # ...)C=O/
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C','c']:  # C(...)C=O/
                aldehyde_number = aldehyde_number + 1

            # aldehyde connected to a ring
            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # C1(...)C=O/
                aldehyde_number = aldehyde_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # C12(...)C=O/
                aldehyde_number = aldehyde_number + 1

            elif s[first_opening_parenthesis-1] == ')':  # ...)(...)C=O/
                index_2 = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index_2)
                if s[second_opening_parenthesis-1] in ['C']:  # C(...)(...)C=O/
                    aldehyde_number = aldehyde_number + 1

        elif len(s) == 3:
            aldehyde_number = aldehyde_number + 1

    # aldehyde as middle characters of the SMILE
    co_index = s.find('(C=O)', 0)
    while co_index != -1:
        if s[co_index-1] in ['C']:  #  C(C=O)...
            aldehyde_number = aldehyde_number + 1

        # aldehyde connected to a ring
        elif s[co_index-2:co_index] in list_tot:  #  C1(C=O)...
            aldehyde_number = aldehyde_number + 1

        elif s[co_index-3:co_index] in list_tot:  #  C12(C=O)...
            aldehyde_number = aldehyde_number + 1

        elif s[co_index-1] == ')':  #  ...)(C=O)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C']:  #  C(...)(C=O)...
                aldehyde_number = aldehyde_number + 1

        co_index = s.find('(C=O)', co_index + 1)

    # aldehyde as last characters of a branch in the SMILE
    co_index = s.find('C=O)', 0)
    while co_index != -1:
        if s[co_index-1] in ['C']:  #  CC=O)...
            aldehyde_number = aldehyde_number + 1

        # aldehyde connected to a ring
        elif s[co_index-2:co_index] in list_tot:  #  C1C=O)...
            aldehyde_number = aldehyde_number + 1

        elif s[co_index-3:co_index] in list_tot:  #  C12C=O)...
            aldehyde_number = aldehyde_number + 1

        elif s[co_index-1] == ')':  #  )C=O)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C']:  #  C(...)C=O)...
                aldehyde_number = aldehyde_number + 1

            # aldehyde connected to a ring
            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)C=O)...
                aldehyde_number = aldehyde_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)C=O)...
                aldehyde_number = aldehyde_number + 1

        co_index = s.find('C=O)', co_index + 1)

    return aldehyde_number

###################   Ketone  #######################################################
def ketone_group(s):
    ketone_number = 0

    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings

    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # ketone as first characters of the SMILE:
    if s.startswith('O=C(C') or s.startswith('O=C(c'):  # /O=C(C... # /O=C(c...
        index = s.find('(', 0)
        first_closing_parenthesis = find_closing_parenthesis(s, index)
        if s[first_closing_parenthesis+1] in ['C','c']:  # /O=C(C...)C or # /O=C(C...)c
            ketone_number = ketone_number + 1

    if s.startswith('O=C1'):  # /O=C1
        if s[4] == 'C': # /O=C1C
            index_1 = s.find('1', 4)
            if s[index_1-1] == 'C':  # /O=C1C...C1
                ketone_number = ketone_number + 1

    # ketone as last character of the SMILE:
    if s[-2:] == '=O':  #  =O/
        # ring-connected
        if s[-4:-2] in list_tot:  #  C1=O/
            if s[-5] == 'C':  #  CC1=O/
                ring_index = s.find(s[-3], 0, -4)
                if s[ring_index-1] in ['C','c']:  #  C1...CC1=O/
                    ketone_number = ketone_number + 1

            elif s[-5] == ')':  #  )C1=O/
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  #  C(...)C1=O/
                    ring_index = s.find(s[-3], 0, first_opening_parenthesis)
                    if s[ring_index-1] in ['C','c']:  #  C1...C(...)C1=O/
                        ketone_number = ketone_number + 1

                elif s[first_opening_parenthesis-1] == ')':  #  )(...)C1=O/
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)C1=O/
                        ring_index = s.find(s[-3], 0, second_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)C1=O/
                            ketone_number = ketone_number + 1

        elif s[-5:-2] in list_tot:  #  C12=O/
            if s[-6] == 'C':  #  CC12=O/
                ring_index = s.find(s[-4:-2], 0, -4)
                if s[ring_index-1] in ['C','c']:  #  C12...CC12=O/
                    ketone_number = ketone_number + 1

            elif s[-6] == ')':  #  )C12=O/
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  #  C(...)C12=O/
                    ring_index = s.find(s[-4:-2], 0, first_opening_parenthesis)
                    if s[ring_index-1] in ['C','c']:  #  C12...C(...)C12=O/
                        ketone_number = ketone_number + 1

                elif s[first_opening_parenthesis-1] == ')':  #  )(...)C12=O/
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)C12=O/
                        ring_index = s.find(s[-4:-2], 0, second_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)C12=O/
                            ketone_number = ketone_number + 1


    if s[-3:] == ')=O':  #  )=O/
        index_2 = s.rfind(')', 0, -2)
        first_opening_parenthesis = find_opening_parenthesis(s, index_2)
        if s[first_opening_parenthesis+1] in ['C','c'] and s[first_opening_parenthesis-1] == 'C':  #  C(C...)=O/  or  C(c...)=O/
            if s[first_opening_parenthesis-2] == 'C':  #  CC(C...)=O/
                ketone_number = ketone_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1C(C...)=O/  or  c1C(C...)=O/
                ketone_number = ketone_number + 1

            elif s[first_opening_parenthesis-4:first_opening_parenthesis-1] in list_tot:  #  C12C(C...)=O/  or  c12C(C...)=O/
                ketone_number = ketone_number + 1

            elif s[first_opening_parenthesis-2] == ')':  #  )C(C...)=O/
                index_3 = s.rfind(')', 0, first_opening_parenthesis-1)
                second_opening_parenthesis = find_opening_parenthesis(s, index_3)
                if s[second_opening_parenthesis-1] == 'C':  #  C(...)C(C...)=O/
                    ketone_number = ketone_number + 1

                elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)C(C...)=O/
                    ketone_number = ketone_number + 1

                elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)C(C...)=O/
                    ketone_number = ketone_number + 1

                elif s[second_opening_parenthesis-1] == ')':  #  )(...)C(C...)=O/
                    index_4 = s.rfind(')', 0, second_opening_parenthesis)
                    third_opening_parenthesis = find_opening_parenthesis(s, index_4)
                    if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)C(C...)=O/
                        ketone_number = ketone_number + 1


    # ketone as middle character of the SMILE:
    co_index = s.find('C(=O)', 0)  #  C(=O)
    while co_index != -1:
        if s[co_index+5] in ['C','c']:  #  C(=O)C  or  C(=O)c
            if s[co_index-1] == 'C' and co_index-1 != -1:  #  CC(=O)C
                ketone_number = ketone_number + 1

            elif s[co_index-2:co_index] in list_tot:  #  C1C(=O)C
                ketone_number = ketone_number + 1

            elif s[co_index-3:co_index] in list_tot:  #  C12C(=O)C
                ketone_number = ketone_number + 1

            elif s[co_index-1] == '(':  #  (C(=O)C
                if s[co_index-2] == 'C':  #  C(C(=O)C
                    ketone_number = ketone_number + 1

                elif s[co_index-3:co_index-1] in list_tot:  #  C1(C(=O)C
                    ketone_number = ketone_number + 1

                elif s[co_index-4:co_index-1] in list_tot:  #  C12(C(=O)C
                    ketone_number = ketone_number + 1

                elif s[co_index-2] == ')':  #  )(C(=O)C
                    index_5 = s.rfind(')', 0, co_index-1)
                    forth_opening_parenthesis = find_opening_parenthesis(s, index_5)
                    if s[forth_opening_parenthesis-1] == 'C':  #  C(...)(C(=O)C
                        ketone_number = ketone_number + 1

            elif s[co_index-1] == ')':  #  )C(=O)C
                index_6 = s.rfind(')', 0, co_index)
                fifth_opening_parenthesis = find_opening_parenthesis(s, index_6)
                if s[fifth_opening_parenthesis-1] == 'C':  #  C(...)C(=O)C
                    ketone_number = ketone_number + 1

                elif s[fifth_opening_parenthesis-2:fifth_opening_parenthesis] in list_tot:  #  C1(...)C(=O)C
                    ketone_number = ketone_number + 1

                elif s[fifth_opening_parenthesis-3:fifth_opening_parenthesis] in list_tot:  #  C12(...)C(=O)C
                    ketone_number = ketone_number + 1

                elif s[fifth_opening_parenthesis-1] == ')':  #  )(...)C(=O)C
                    index_7 = s.rfind(')', 0, fifth_opening_parenthesis)
                    sixth_opening_parenthesis = find_opening_parenthesis(s, index_7)
                    if s[sixth_opening_parenthesis-1] == 'C':  #  C(...)(...)C(=O)C
                        ketone_number = ketone_number + 1

        co_index = s.find('C(=O)C', co_index + 1)


    co_index2 = s.find(')=O)', 0)  #  )=O)
    while co_index2 != -1:
        index_8 = s.rfind(')', 0, co_index2+1)
        seventh_opening_parenthesis = find_opening_parenthesis(s, index_8)
        if s[seventh_opening_parenthesis+1] in ['C','c'] and s[seventh_opening_parenthesis-1] == 'C':  #  C(C...)=O)  or  #  C(c...)=O)
            if s[seventh_opening_parenthesis-2] == 'C': #  CC(C...)=O)  keton is the last characters of a branch
                ketone_number = ketone_number + 1

            elif s[seventh_opening_parenthesis-2] == '(':  #  (C(C...)=O)
                if s[seventh_opening_parenthesis-3] == 'C':  #  C(C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[seventh_opening_parenthesis-4:seventh_opening_parenthesis-2] in list_tot:  #  C1(C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[seventh_opening_parenthesis-5:seventh_opening_parenthesis-2] in list_tot:  #  C12(C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[seventh_opening_parenthesis-3] == ')':  #  )(C(C...)=O)
                    index_9 = s.rfind(')', 0, seventh_opening_parenthesis-2)
                    eighth_opening_parenthesis = find_opening_parenthesis(s, index_9)
                    if s[eighth_opening_parenthesis-1] == 'C':  #  C(...)(C(C...)=O)
                        ketone_number = ketone_number + 1

            # ketone as last characters of a branch
            elif s[seventh_opening_parenthesis-2] == ')':  #  )C(C...)=O)
                index_10 = s.rfind(')', 0, seventh_opening_parenthesis-1)
                ninth_opening_parenthesis = find_opening_parenthesis(s, index_10)
                if s[ninth_opening_parenthesis-1] == 'C':  #  C(...)C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[ninth_opening_parenthesis-2:ninth_opening_parenthesis] in list_tot:  #  C1(...)C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[ninth_opening_parenthesis-3:ninth_opening_parenthesis] in list_tot:  #  C12(...)C(C...)=O)
                    ketone_number = ketone_number + 1

                elif s[ninth_opening_parenthesis-1] == ')':  #  )(...)C(C...)=O)
                    index_11 = s.rfind(')', 0, ninth_opening_parenthesis)
                    tenth_opening_parenthesis = find_opening_parenthesis(s, index_11)
                    if s[tenth_opening_parenthesis-1] == 'C':  #  C(...)(...)C(C...)=O)
                            ketone_number = ketone_number + 1

        co_index2 = s.find(')=O)', co_index2 + 1)


    co_index2 = s.find('=O)', 0)  #  =O)
    while co_index2 != -1:
        # ring-connected
        if s[co_index2-2:co_index2] in list_tot:  #  C1=O)
            if s[co_index2-3] == 'C':  #  CC1=O)
                ring_index = s.find(s[co_index2-1], 0, co_index2-2)
                if s[ring_index-1] in ['C','c']:  #  C1...CC1=O)
                    ketone_number = ketone_number + 1

            elif s[co_index2-3] == ')':  #  )C1=O)
                index = s.rfind(')', 0, co_index2-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  #  C(...)C1=O)
                    ring_index = s.find(s[co_index2-1], 0, first_opening_parenthesis)
                    if s[ring_index-1] in ['C','c']:  #  C1...C(...)C1=O)
                        ketone_number = ketone_number + 1

                elif s[first_opening_parenthesis-1] == ')':  #  )(...)C1=O)
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)C1=O)
                        ring_index = s.find(s[co_index2-1], 0, second_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)C1=O)
                            ketone_number = ketone_number + 1


        if s[co_index2-3:co_index2] in list_tot:  #  C12=O)
            if s[co_index2-4] == 'C':  #  CC12=O)
                ring_index = s.find(s[co_index2-2:co_index2], 0, co_index2-2)
                if s[ring_index-1] in ['C','c']:  #  C12...CC12=O)
                    ketone_number = ketone_number + 1

            elif s[co_index2-4] == ')':  #  )C12=O)
                index = s.rfind(')', 0, co_index2-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  #  C(...)C12=O)
                    ring_index = s.find(s[co_index2-2:co_index2], 0, first_opening_parenthesis)
                    if s[ring_index-1] in ['C','c']:  #  C12...C(...)C12=O)
                        ketone_number = ketone_number + 1

                elif s[first_opening_parenthesis-1] == ')':  #  )(...)C12=O)
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)C12=O)
                        ring_index = s.find(s[co_index2-2:co_index2], 0, second_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)C12=O)
                            ketone_number = ketone_number + 1

        co_index2 = s.find('=O)', co_index2 + 1)



    # Alkoxy Radical ketone as first characters of the SMILE:
    allowed_prefixes = ['O=C([O])', '[O]C(=O)']
    for prefix in allowed_prefixes:
        if s.startswith(prefix):  #  /O=C([O])... or /[O]C(=O)...
            ketone_number = ketone_number + 1

    # Alkoxy Radical ketone as last characters of the SMILE:
    if s[-8:] == 'C(=O)[O]' or s[-8:] == 'C([O])=O':  #  ...C(=O)[O]/ or ...C([O])=O/
        ketone_number = ketone_number + 1

    # Alkoxy Radical ketone as middle characters of the SMILE:
    co_index = s.find('C(=O)[O])', 0)
    while co_index != -1:
        ketone_number = ketone_number + 1

        co_index = s.find('C(=O)[O])', co_index + 1)

    co_index = s.find('C([O])=O)', 0)
    while co_index != -1:
        ketone_number = ketone_number + 1

        co_index = s.find('C([O])=O)', co_index + 1)

    if s == 'O=C[O]' or s == '[O]C=O':  # O=C[O]
        ketone_number = ketone_number + 1


    # Peroxy Radical ketone as first characters of the SMILE:
    allowed_prefixes = ['O=C(O[O])', '[O]OC(=O)']
    for prefix in allowed_prefixes:
        if s.startswith(prefix):  #  /O=C(O[O])... or /[O]OC(=O)...
            ketone_number = ketone_number + 1

    # Peroxy Radical ketone as last characters of the SMILE:
    if s[-9:] == 'C(=O)O[O]' or s[-9:] == 'C(O[O])=O':  #  ...C(=O)O[O]/ or ...C(O[O])=O/
        ketone_number = ketone_number + 1

    # Peroxy Radical ketone as middle characters of the SMILE:
    co_index = s.find('C(=O)O[O])', 0)
    while co_index != -1:
        ketone_number = ketone_number + 1

        co_index = s.find('C(=O)O[O])', co_index + 1)

    co_index = s.find('C(O[O])=O)', 0)
    while co_index != -1:
        ketone_number = ketone_number + 1

        co_index = s.find('C(O[O])=O)', co_index + 1)

    if s == 'O=CO[O]' or s == '[O]OC=O':  # O=CO[O]
        ketone_number = ketone_number + 1



    return ketone_number

###################   Carboxilic acid  #######################################################
def carboxylic_acid_group(s):
    carboxylic_acid_number = 0
    # carboxylic_acid as first characters of the SMILES:
    allowed_prefixes = ['O=C(O)', 'OC(=O)']
    for prefix in allowed_prefixes:
        if s.startswith(prefix):  #  /O=C(O)... or /OC(=O)...
            carboxylic_acid_number = carboxylic_acid_number + 1

    # carboxylic_acid as last characters of the SMILES:
    if s[-6:] == 'C(=O)O' or s[-6:] == 'C(O)=O':  #  ...C(=O)O/ or ...C(O)=O/
        carboxylic_acid_number = carboxylic_acid_number + 1

    # carboxylic_acid as middle characters of the SMILES:
    co_index = s.find('C(=O)O)', 0)
    while co_index != -1:
        carboxylic_acid_number = carboxylic_acid_number + 1

        co_index = s.find('C(=O)O)', co_index + 1)

    co_index = s.find('C(O)=O)', 0)
    while co_index != -1:
        carboxylic_acid_number = carboxylic_acid_number + 1

        co_index = s.find('C(O)=O)', co_index + 1)

    if s == 'O=CO' or s == 'OC=O':
        carboxylic_acid_number = carboxylic_acid_number + 1

    return carboxylic_acid_number

###################   Ester  #######################################################
def ester_group(s):
    ester_number = 0
    nitroester_number = 0
    cyc_number = find_highest_digit(s)

    digits = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

    # Replace '[N+](=O)[O-]', 'O=[N+][O-]', and '[O-][N+](=O)' with 'N(=O)=O' and 'O=N(=O)' respectively:
    s = s.replace('[N+](=O)[O-]', 'N(=O)=O')
    s = s.replace('O=[N+][O-]', 'O=N(=O)')
    s = s.replace('[O-][N+](=O)', 'O=N(=O)')

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings
    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # ester as first characters of the SMILE
    if s[0:5] == 'O=C(O':  # /O=C(O
        if s[5] in ['C','c']:  # /O=C(OC or # /O=C(Oc
            index = s.find('(', 3)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1] in ['C','c']:  # /O=C(OC...)C
                co_index = s.find('N(=O)=O', first_closing_parenthesis+1)
                if co_index != -1:  # /O=C(OC...)C...N(=O)=O...
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

    elif s[0:4] == 'O=C(':
        if s[4] in ['C','c']:  #  /O=C(C or /O=C(c
            index = s.find('(', 3)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1:first_closing_parenthesis+3] in ['OC','Oc']:  #  /O=C(C...)OC
                co_index = s.find('N(=O)=O', 3, first_closing_parenthesis)
                if co_index != -1:  #  /O=C(C...N(=O)=O...)OC
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1


    elif s[0:4] == 'O=C1':
        if s[4] in ['C','c']:  #  /O=C1C or /O=C1c
            index = s.find('1', 4)
            if s[index-1] == 'O':  #  /O=C1C...O1
                ester_number = ester_number + 1

            elif s[index-1] in digits:
                if s[index-2] == 'O':  #  /O=C1C...O21
                    ester_number = ester_number + 1

        elif s[4] == 'O':  #  /O=C1O or /O=C1O
            index = s.find('1', 4)
            if s[index-1] in ['C','c']:  #  /O=C1O...C1
                ester_number = ester_number + 1

            elif s[index-1] in digits:
                if s[index-2] == 'O':  #  /O=C1O...C21
                    ester_number = ester_number + 1


    # Ester as last characters of the SMILE
    if s[-2:] == '=O':  #  =O/
        if s[-3] == ')':    #  )=O/
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-2:first_opening_parenthesis+2] in ['OC(C', 'OC(c']:  #  OC(C...)=O/
                if s[first_opening_parenthesis-3] in ['C','c'] :  #  COC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  COC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1OC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  C1OC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12OC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  C12OC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )OC(C...)=O
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(C...)=O
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)OC(C...)=O/
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C1(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)OC(C...)=O/
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C1(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1


                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)OC(C...)=O/
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)OC(C...)=O/
                            co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                            if co_index != -1:
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

            elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] in ['C(OC','C(Oc']:  #  C(OC...)=O/
                if s[first_opening_parenthesis-2] in ['C','c']:  #  CC(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...CC(OC...)=O/ or O=N(=O)...CC(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1C(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...C1C(OC...)=O/ or O=N(=O)...C1C(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-1] in list_tot:  #  C12C(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...C12C(OC...)=O/ or O=N(=O)...C12C(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2] == ')':  #  )C(OC...)=O/
                    index = s.rfind(')', 0, first_opening_parenthesis-1)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C1(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C1(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C12(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C12(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)C(OC...)=O/
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)C(OC...)=O/
                            co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                            co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                            if co_index != -1 or co_index1 != -1:  #  C(...)(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C(...)(...)C(OC...)=O/
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1


        elif s[-4:-2] == 'OC':  # ...OC=O
            if len(s) >= 5 and s[-5] == 'C':  # ...COC=O
                ester_number = ester_number + 1

            elif s[-6:-4] in list_tot:  # ...C1OC=O  or  # ...c1OC=O
                ester_number = ester_number + 1

            elif s[-7:-4] in list_tot:  # ...C12OC=O  or  # ...c12OC=O
                ester_number = ester_number + 1

            elif len(s) >= 5 and s[-5] == ')':  # ...)OC=O
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  # ...C(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # ...C1(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # ...C12(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OC=O
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)OC=O
                        ester_number = ester_number + 1


        elif s[-4:-2] in list_tot:  # ...C1=O/
            if s[-5] == 'O': # ...OC1=O/
                if s[-6] == 'C':  # ...COC1=O/
                    ring_index = s.find(s[-3], 0, -5)
                    if s[ring_index-1] in ['C','c']:  #  C1...COC1=O/
                        ester_number = ester_number + 1

                elif s[-6] == ')':  # ...)OC1=O
                    index = s.rfind(')', 0, -5)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC1=O/
                        ring_index = s.find(s[-3], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)OC1=O/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC1=O/
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC1=O/
                            ring_index = s.find(s[-3], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)OC1=O/
                                ester_number = ester_number + 1



        elif s[-5:-2] in list_tot:  # ...C12=O/
            if s[-6] == 'O': # ...OC12=O/
                if s[-7] == 'C':  # ...COC21=O/
                    ring_index = s.find(s[-4:-2], 0, -5)
                    if s[ring_index-1] in ['C','c']:  #  C12...COC12=O/
                        ester_number = ester_number + 1

                elif s[-7] == ')':  # ...)OC12=O
                    index = s.rfind(')', 0, -5)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC12=O/
                        ring_index = s.find(s[-4:-2], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)OC1=O/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC12=O/
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC12=O/
                            ring_index = s.find(s[-4:-2], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)OC12=O/
                                ester_number = ester_number + 1



    # Ester as last characters of a channel in the SMILE
    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-1] == ')':  # ...)=O)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-2:first_opening_parenthesis+2] in ['OC(C', 'OC(c']:  #  OC(C...)=O)
                if s[first_opening_parenthesis-3] in ['C','c'] :  #  COC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  COC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1OC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  C1OC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12OC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  C12OC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )OC(C...)=O)
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C1(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C1(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)OC(C...)=O)
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)OC(C...)=O)
                            co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                            if co_index1 != -1:  #  C(...)(...)OC(C...N(=O)=O...)=O)
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

            elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] in ['C(OC','C(Oc']:  #  C(OC...)=O)
                if s[first_opening_parenthesis-2] in ['C','c']:  #  CC(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...CC(OC...)=O) or ...N(=O)=O...CC(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1C(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C1C(OC...)=O) or ...N(=O)=O...C1C(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-1] in list_tot:  #  C12C(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C12C(OC...)=O) or ...N(=O)=O...C12C(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2] == ')':  #  )C(OC...)=O)
                    index = s.rfind(')', 0, first_opening_parenthesis-1)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C(...)C(OC...)=O) or ...N(=O)=O...C(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C1(...)C(OC...)=O) or ...N(=O)=O...C1(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1


                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C12(...)C(OC...)=O) or ...N(=O)=O...C1(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)C(OC...)=O)
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)C(OC...)=O)
                            co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                            co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                            if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C(...)(...)C(OC...)=O) or ...N(=O)=O...C(...)(...)C(OC...)=O)
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

        elif s[co_index-2:co_index] == 'OC':  # ...OC=O)
            if s[co_index-3] == 'C':  # ...COC=O)
                ester_number = ester_number + 1

            elif s[co_index-4:co_index-2] in list_tot:  # ...C1OC=O)  or  # ...c1OC=O)
                ester_number = ester_number + 1

            elif s[co_index-5:co_index-2] in list_tot:  # ...C12OC=O)  or  # ...c12OC=O)
                ester_number = ester_number + 1

            elif s[co_index-3] == ')':  # ...)OC=O)
                index = s.rfind(')', 0, co_index-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  # ...C(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # ...C1(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # ...C12(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OC=O)
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)OC=O)
                        ester_number = ester_number + 1


        elif s[co_index-2:co_index] in list_tot:  # ...C1=O)
            if s[co_index-3] == 'O': # ...OC1=O)
                if s[co_index-4] == 'C':  # ...COC1=O)
                    ring_index = s.find(s[co_index-1], 0, co_index-4)
                    if s[ring_index-1] in ['C','c']:  #  C1...COC1=O)
                        ester_number = ester_number + 1

                elif s[co_index-4] == ')':  # ...)OC1=O)
                    index = s.rfind(')', 0, co_index-3)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC1=O)
                        ring_index = s.find(s[co_index-1], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)OC1=O)
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC1=O)
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC1=O)
                            ring_index = s.find(s[co_index-1], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)OC1=O)
                                ester_number = ester_number + 1


        elif s[co_index-3:co_index] in list_tot:  # ...C12=O)
            if s[co_index-4] == 'O': # ...OC12=O)
                if s[co_index-5] == 'C':  # ...COC12=O)
                    ring_index = s.find(s[co_index-2:co_index], 0, co_index-4)
                    if s[ring_index-1] in ['C','c']:  #  C12...COC12=O)
                        ester_number = ester_number + 1

                elif s[co_index-5] == ')':  # ...)OC12=O)
                    index = s.rfind(')', 0, co_index-3)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC12=O)
                        ring_index = s.find(s[co_index-2:co_index], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)OC12=O)
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC12=O)
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC12=O)
                            ring_index = s.find(s[co_index-2:co_index], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)OC12=O)
                                ester_number = ester_number + 1



        co_index = s.find('=O)', co_index + 1)

    # Ester as middle characters in the SMILE:
    co_index1 = s.find('C(=O)O', 0)
    while co_index1 != -1:
        if len(s) > co_index1+6 and s[co_index1+6] in ['C','c']:  # C(=O)OC
            if s[co_index1-1] in ['C','c'] and co_index1-1 != -1:  # CC(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-2:co_index1] in list_tot and co_index1-2 != -2:  # C1C(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-3:co_index1] in list_tot  and co_index1-3 != -3:  # C12C(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-1] == ')':  # )C(=O)OC
                index = s.rfind(')', 0, co_index1)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] in ['C', 'c']:  # C(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # C1(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # C12(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] in [')']:  # )(...)C(=O)OC
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']: # C(...)(...)C(=O)OC
                        co_index2 = s.find('N(=O)=O', 0, co_index1)
                        co_index3 = s.find('O=N(=O)', 0, co_index1)
                        if co_index2 != -1 or co_index3 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

            elif s[co_index1-1] == '(':  # (C(=O)OC
                if s[co_index1-2] in ['C','c']:  # C(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-3:co_index1-1] in list_tot:  # C1(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-4:co_index1-1] in list_tot:  # C12(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-2] == ')':  # )(C(=O)OC
                    index = s.rfind(')', 0, co_index1-1)
                    third_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[third_opening_parenthesis-1] == 'C':  # C(...)(C(=O)OC
                        index2 = s.find('(', co_index1-1)
                        first_closing_parenthesis = find_closing_parenthesis(s, index2)
                        co_index2 = s.find('N(=O)=O', 0, co_index1)
                        co_index3 = s.find('O=N(=O)', 0, co_index1)
                        co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                        if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1


        co_index1 = s.find('C(=O)O', co_index1 + 1)

    co_index2 = s.find(')=O)', 0)
    while co_index2 != -1:
        index = s.rfind(')', 0, co_index2+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-2:first_opening_parenthesis+3] in ['(C(OC','(C(Oc']:  # (C(OC...)=O)
            if s[first_opening_parenthesis-3] in ['C','c']:  # C(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  # C1(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  # C12(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-3] == ')': # )(C(OC...)=O)
                index = s.rfind(')', 0, first_opening_parenthesis-2)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:   # C(...)(C(OC...)=O)
                    co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    co_index5 = s.find('N(=O)=O', co_index2)
                    if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

        elif s[first_opening_parenthesis-3:first_opening_parenthesis+2] in ['(OC(C','(OC(c']:
            if s[first_opening_parenthesis-4] in ['C','c']:  #  C(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-5:first_opening_parenthesis-3] in list_tot:  #  C1(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-6:first_opening_parenthesis-3] in list_tot:  #  C12(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-4] == ')':  #  )(OC(C...)=O)
                index = s.rfind(')', 0, first_opening_parenthesis-3)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)(OC(C...)=O)
                    co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                    if co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

        co_index2 = s.find(')=O)', co_index2 + 1)

    co_index3 = s.find('OC(=O)', 0)
    if co_index3 != 0:  #  OC(=O)C...
        while co_index3 != -1:
            if s[co_index3+6] in ['C','c']:  #  OC(=O)C...
                if s[co_index3-1] in ['C','c'] and co_index3-1 != -1:  #  COC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-2:co_index3] in list_tot and co_index3-2 != -2:  #  C1OC(=O)C...  not to be confused with OC(=O)C...C1/
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-3:co_index3] in list_tot and co_index3-3 != -3:  #  C12OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-1] == ')':  #  )OC(=O)C...
                    index = s.rfind(')', 0, co_index3)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] in [')']:  #  )(...)OC(=O)C...
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] in ['C']:  #  C(...)(...)OC(=O)C...
                            co_index4 = s.find('N(=O)=O', co_index3)
                            if co_index4 != -1:
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

                elif s[co_index3-1] == '(':  #  (OC(=O)C...
                    index = s.find('(', co_index3-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[co_index3-2] in ['C','c']:  #  C(OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[co_index3-3:co_index3-1] in list_tot:  #  C1(OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[co_index3-4:co_index3-1] in list_tot:  #  C12(OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[co_index3-2] == ')':  #  )(OC(=O)C...
                        index = s.rfind(')', 0, co_index3-1)
                        first_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[first_opening_parenthesis-1] == 'C':  #  C(...)(OC(=O)C...
                            co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                            if co_index4 != -1:
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

            co_index3 = s.find('OC(=O)', co_index3 + 1)


    if cyc_number != None: # if there is any cycle
    # Ending with cycle
        for i in range(1, cyc_number+1):
            co_index = s.find('(OC' + str(i) + '=O)', 0)
            while co_index != -1:  # ...(OC1=O)...
                if s[co_index-1] == 'C':  # ...C(OC1=O)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...C(OC1=O)...
                        ester_number = ester_number + 1

                elif s[co_index-1] == ')':  # ...)(OC1=O)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)(OC1=O)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(OC1=O)...
                            ester_number = ester_number + 1

                co_index = s.find('(OC' + str(i) + '=O)', co_index + 1)


        for i in range(1, cyc_number+1):
            co_index = s.find('(C(=O)O' + str(i) + ')', 0)
            while co_index != -1:  # ...(C(=O)O1)...
                if s[co_index-1] == 'C':  # ...C(C(=O)O1)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...C(C(=O)O1)...
                        ester_number = ester_number + 1

                elif s[co_index-1] == ')':  # ...)(C(=O)O1)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)(C(=O)O1)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(C(=O)O1)...
                            ester_number = ester_number + 1

                co_index = s.find('(C(=O)O' + str(i) + ')', co_index + 1)


        for i in range(1, cyc_number+1):
            co_index = s.find('(C(O' + str(i) + ')=O)', 0)
            while co_index != -1:  # ...(C(O1)=O)...
                if s[co_index-1] == 'C':  # ...C(C(O1)=O)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...C(C(O1)=O)...
                        ester_number = ester_number + 1

                elif s[co_index-1] == ')':  # ...)(C(O1)=O)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)(C(O1)=O)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(C(O1)=O)...
                            ester_number = ester_number + 1

                co_index = s.find('(C(O' + str(i) + ')=O)', co_index + 1)


        for i in range(1, cyc_number+1):
            co_index = s.find('C(=O)O' + str(i) + ')', 0)
            while co_index != -1:  # ...C(=O)O1)...
                if s[co_index-1] == 'C':  # ...CC(=O)O1)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...CC(=O)O1)...
                        ester_number = ester_number + 1

                elif s[co_index-1] == ')':  # ...)C(=O)O1)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)C(=O)O1)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)C(=O)O1)...
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)C(=O)O1)...
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)C(=O)O1)...
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)C(=O)O1)...
                                ester_number = ester_number + 1

                co_index = s.find('(C(=O)O' + str(i) + ')', co_index + 1)


        for i in range(1, cyc_number+1):
            co_index = s.find('C(O' + str(i) + ')=O)', 0)
            while co_index != -1:  # ...C(O1)=O)...
                if s[co_index-1] == 'C':  # ...CC(O1)=O)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...CC(O1)=O)...
                        ester_number = ester_number + 1

                elif s[co_index-1] == ')':  # ...)C(O1)=O)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)C(O1)=O)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)C(O1)=O)...
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)C(O1)=O)...
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)C(O1)=O)...
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)C(O1)=O)...
                                ester_number = ester_number + 1

                co_index = s.find('C(O' + str(i) + ')=O)', co_index + 1)


        for i in range(1, cyc_number+1):
            end_str = 'C(=O)O' + str(i)
            if s.endswith(end_str):
                if len(s) > len(end_str) and s[-len(end_str) - 1] == 'C':   # CC(=O)O1/
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...CC(=O)O1/
                        ester_number = ester_number + 1

                elif len(s) > len(end_str) and s[-len(end_str) - 1] == ')':   # )C(=O)O1/
                    index1 = s.rfind(')', 0, -len(end_str))
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)C(=O)O1/
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)C(=O)O1/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)C(=O)O1/
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)C(=O)O1/
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)C(=O)O1/
                                ester_number = ester_number + 1


        for i in range(1, cyc_number+1):
            end_str = 'C(O' + str(i) + ')=O'
            if s.endswith(end_str):
                if len(s) > len(end_str) and s[-len(end_str) - 1] == 'C':   # CC(O1)=O/
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...CC(O1)=O/
                        ester_number = ester_number + 1

                elif len(s) > len(end_str) and s[-len(end_str) - 1] == ')':   # )C(O1)=O/
                    index1 = s.rfind(')', 0, -len(end_str))
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)C(O1)=O/
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)C(O1)=O/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)C(O1)=O/
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)C(O1)=O/
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)C(O1)=O/
                                ester_number = ester_number + 1


    return ester_number

###################   Ether  #######################################################
def ether_group(s):
    ether_number = 0
    cyc_number = find_highest_digit(s)
    co_index = s.find('OC', 0)
    while co_index != -1:
        if s[co_index-1] in ['C'] and co_index != 0:  # COC   avoin /OC...C to be counted as ether
            # excluding alicyclic ether:
            if cyc_number != None:  # if there is any ring
                i = 1
                open_loop = True
                while i <= cyc_number and open_loop:
                    A = s.find(str(i), 0)
                    B = s.find(str(i), A+1)
                    if A < co_index-1 and B > co_index+1:    # 1...COC...1
                        D = s.find('(', A, B)
                        if D == -1:  # ether is in the ring
                            break

                        else:
                            while D != -1:
                                E = find_closing_parenthesis(s, D)
                                if D < co_index-1 and co_index+1 < E < B:  # 1...(...COC...)...1  ether is outside of the ring
                                    # excluding ester:
                                    if s[co_index+2:co_index+4] != '=O' and s[co_index-3:co_index-1] != 'O=':  #  ...COC=O/ #  /O=COC...
                                        if len(s) > co_index+2 and s[co_index+2] == '(': #  ...COC(...
                                            index = s.find('(', co_index+2)
                                            first_closing_parenthesis = find_closing_parenthesis(s, index)
                                            if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...COC(=O)=O
                                                ether_number = ether_number + 1
                                                open_loop = False
                                                break

                                            else:
                                                D = s.find('(', D + 1, B)
                                                if D == -1:  # No more opening parenthesis and ether is in the ring
                                                    open_loop = False
                                                    break

                                        else:
                                            ether_number = ether_number + 1
                                            open_loop = False
                                            break

                                    else:
                                        D = s.find('(', D + 1, B)
                                        if D == -1:  # No more opening parenthesis and ether is in the ring
                                            open_loop = False
                                            break

                                else:
                                    D = s.find('(', D + 1, B)
                                    if D == -1:  # No more opening parenthesis and ether is in the ring
                                        open_loop = False
                                        break

                    else:
                        if i == cyc_number:  # last ring
                            # excluding ester:
                            if s[co_index+2:co_index+4] != '=O' and s[co_index-3:co_index-1] != 'O=':  #  ...COC=O/ #  /O=COC...
                                if len(s) > co_index+2 and s[co_index+2] == '(': #  ...COC(...
                                    index = s.find('(', co_index+2)
                                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...COC(=O)
                                        ether_number = ether_number + 1
                                        open_loop = False
                                        break

                                    else:
                                        i = i + 1

                                else:
                                    ether_number = ether_number + 1
                                    open_loop = False
                                    break

                            else:
                                i = i + 1

                        else:
                            i = i + 1

            else:
                # excluding ester:
                if s[co_index+2:co_index+4] != '=O' and s[co_index-3:co_index-1] != 'O=':  #  ...COC=O/ #  /O=COC...
                    if len(s) > co_index+2 and s[co_index+2] == '(': #  ...COC(...
                        index = s.find('(', co_index+2)
                        first_closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...COC(=O)
                            ether_number = ether_number + 1

                    else:
                        ether_number = ether_number + 1


        elif s[co_index-1] in ['(']:  # ...(OC...
            if s[co_index-2] in ['C']:  #  ...C(OC...
                # excluding alicyclic ether:
                if cyc_number != None:  # if there is any ring
                    i = 1
                    open_loop = True
                    while i <= cyc_number and open_loop:
                        A = s.find(str(i), 0)
                        B = s.find(str(i), A+1)
                        C = find_closing_parenthesis(s, co_index-1)
                        if A < co_index-2 and co_index-1 < B < C:  # 1...C(OC...1...)  ether is in the ring
                            open_loop = False
                            break

                        else:
                            if i == cyc_number:  # last ring
                                # excluding ester:
                                index = s.find('(', co_index-1)
                                first_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+2:co_index+4] != '=O' and s[co_index-4:co_index-2] != 'O=' and s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  ...C(OC=O) #  /O=C(OC...)=O
                                    if s[co_index+2] == '(': #  ...C(OC(...
                                        index = s.find('(', co_index+2)
                                        second_closing_parenthesis = find_closing_parenthesis(s, index)
                                        if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)=O)=O
                                            ether_number = ether_number + 1
                                            open_loop = False
                                            break

                                        else:
                                            i = i + 1

                                    else:
                                        ether_number = ether_number + 1
                                        open_loop = False
                                        break

                                else:
                                    i = i + 1

                            else:
                                i = i + 1

                else:
                    # excluding ester:
                    index = s.find('(', co_index-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[co_index+2:co_index+4] != '=O' and s[co_index-4:co_index-2] != 'O=' and s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  ...C(OC=O) #  /O=C(OC...)=O
                        if s[co_index+2] == '(': #  ...C(OC(...
                            index = s.find('(', co_index+2)
                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                            if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                ether_number = ether_number + 1

                        else:
                            ether_number = ether_number + 1


            elif s[co_index-2] in [')']:  #  ...)(OC...
                index = s.rfind(')', 0, co_index-1)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':    #C(...)(OC...
                    # excluding alicyclic ether:
                    if cyc_number != None:  # if there is any ring
                        i = 1
                        open_loop = True
                        while i <= cyc_number and open_loop:
                            A = s.find(str(i), 0)
                            B = s.find(str(i), A+1)
                            C = find_closing_parenthesis(s, co_index-1)
                            if A < first_opening_parenthesis-1 and co_index-1 < B < C:  # 1...C(...)(OC...1...)  ether is in the ring
                                open_loop = False
                                break

                            else:
                                if i == cyc_number:  # last ring
                                    # excluding ester:
                                    index = s.find('(', co_index-1)
                                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[co_index+2:co_index+4] != '=O' and s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O' and s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  ...C(OC=O) #  /O=C(OC...)=O
                                        if s[co_index+2] == '(': #  ...C(OC(...
                                            index = s.find('(', co_index+2)
                                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                                            if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                                ether_number = ether_number + 1
                                                open_loop = False
                                                break

                                            else:
                                                i = i + 1

                                        else:
                                            ether_number = ether_number + 1
                                            open_loop = False
                                            break

                                    else:
                                        i = i + 1

                                else:
                                    i = i + 1

                    else:
                        # excluding ester:
                        index = s.find('(', co_index-1)
                        first_closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[co_index+2:co_index+4] != '=O' and s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O' and s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  ...C(OC=O) #  /O=C(OC...)=O
                            if s[co_index+2] == '(': #  ...C(OC(...
                                index = s.find('(', co_index+2)
                                second_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                    ether_number = ether_number + 1

                            else:
                                ether_number = ether_number + 1


        elif s[co_index-1] in [')']:  #  ...)OC...
            index2 = s.rfind(')', 0, co_index)
            second_opening_parenthesis = find_opening_parenthesis(s, index2)
            if s[second_opening_parenthesis-1] == 'C':  #  ...C(...)OC...
                # excluding alicyclic ether:
                if cyc_number != None:  # if there is any ring
                    i = 1
                    open_loop = True
                    while i <= cyc_number and open_loop:
                        A = s.find(str(i), 0)
                        B = s.find(str(i), A+1)
                        if A < second_opening_parenthesis-1 and B > co_index+1:    # 1...C(...)OC...1
                            D = s.find('(', A, B)
                            if D == -1:  # ether is in the ring
                                break

                            else:
                                while D != -1:
                                    E = find_closing_parenthesis(s, D)
                                    if D < second_opening_parenthesis-1 and co_index+1 < E < B:  # 1...(...C(...)OC...)...1  ether is out of the ring
                                        # excluding ester:
                                        if s[co_index+2:co_index+4] != '=O' and s[second_opening_parenthesis+1:second_opening_parenthesis+3] != '=O' and s[second_opening_parenthesis-3:second_opening_parenthesis-1] != 'O=':  #  ...C(OC=O) #  /O=C(OC...)=O
                                            if len(s) > co_index+2 and s[co_index+2] == '(': #  ...C)OC(...
                                                index = s.find('(', co_index+2)
                                                second_closing_parenthesis = find_closing_parenthesis(s, index)
                                                if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                                    ether_number = ether_number + 1
                                                    open_loop = False
                                                    break

                                                else:
                                                    D = s.find('(', D + 1, B)
                                                    if D == -1:  # No more opening parenthesis and ether is in the ring
                                                        open_loop = False
                                                        break

                                            else:
                                                ether_number = ether_number + 1
                                                open_loop = False
                                                break

                                        else:
                                            D = s.find('(', D + 1, B)
                                            if D == -1:  # No more opening parenthesis and ether is in the ring
                                                open_loop = False
                                                break

                                    else:
                                        D = s.find('(', D + 1, B)
                                        if D == -1:  # No more opening parenthesis and ether is in the ring
                                            open_loop = False
                                            break

                        else:
                            if i == cyc_number:  # last ring
                                # excluding ester:
                                if s[co_index+2:co_index+4] != '=O' and s[second_opening_parenthesis+1:second_opening_parenthesis+3] != '=O' and s[second_opening_parenthesis-3:second_opening_parenthesis-1] != 'O=':  #  ...C(OC=O) #  /O=C(OC...)=O
                                    if len(s) > co_index+2 and s[co_index+2] == '(': #  ...C)OC(...
                                        index = s.find('(', co_index+2)
                                        second_closing_parenthesis = find_closing_parenthesis(s, index)
                                        if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                            ether_number = ether_number + 1
                                            open_loop = False
                                            break

                                        else:
                                            i = i + 1

                                    else:
                                        ether_number = ether_number + 1
                                        open_loop = False
                                        break

                                else:
                                    i = i + 1

                            else:
                                i = i + 1

                else:
                    # excluding ester:
                    if s[co_index+2:co_index+4] != '=O' and s[second_opening_parenthesis+1:second_opening_parenthesis+3] != '=O' and s[second_opening_parenthesis-3:second_opening_parenthesis-1] != 'O=':  #  ...C(OC=O) #  /O=C(OC...)=O
                        if len(s) > co_index+2 and s[co_index+2] == '(': #  ...C)OC(...
                            index = s.find('(', co_index+2)
                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                            if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...C(OC(=O)
                                ether_number = ether_number + 1

                        else:
                            ether_number = ether_number + 1


            elif s[second_opening_parenthesis-1] == ')':  #  ...)(...)OC...
                index3 = s.rfind(')', 0, second_opening_parenthesis)
                third_opening_parenthesis = find_opening_parenthesis(s, index3)
                if s[third_opening_parenthesis-1] == 'C':   #  ...C(...)(...)OC...
                    # excluding alicyclic ether:
                    if cyc_number != None:  # if there is any ring
                        i = 1
                        open_loop = True
                        while i <= cyc_number and open_loop:
                            A = s.find(str(i), 0)
                            B = s.find(str(i), A+1)
                            if A < third_opening_parenthesis-1 and B > co_index+1:    # 1...C(...)(...)OC...1
                                D = s.find('(', A, B)
                                if D == -1:  # ether is in the ring
                                    break

                                else:
                                    while D != -1:
                                        E = find_closing_parenthesis(s, D)
                                        if D < third_opening_parenthesis-1 and co_index+1 < E < B:  # 1...(...C(...)(...)OC...)...1  ether is out of the ring
                                            # excluding ester:
                                            if s[co_index+2:co_index+4] != '=O':
                                                if len(s) > co_index+2 and s[co_index+2] == '(': #  C(...)(...)OC(...
                                                    index = s.find('(', co_index+2)
                                                    second_closing_parenthesis = find_closing_parenthesis(s, index)
                                                    if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+2:co_index+4] != '=O':
                                                        ether_number = ether_number + 1
                                                        open_loop = False
                                                        break

                                                    else:
                                                        D = s.find('(', D + 1, B)
                                                        if D == -1:  # No more opening parenthesis and ether is in the ring
                                                            open_loop = False
                                                            break

                                                else:
                                                    ether_number = ether_number + 1
                                                    open_loop = False
                                                    break

                                            else:
                                                D = s.find('(', D + 1, B)
                                                if D == -1:  # No more opening parenthesis and ether is in the ring
                                                    open_loop = False
                                                    break


                                        else:
                                            D = s.find('(', D + 1, B)
                                            if D == -1:  # No more opening parenthesis and ether is in the ring
                                                open_loop = False
                                                break


                            else:
                                if i == cyc_number:  # last ring
                                    # excluding ester:
                                    if s[co_index+2:co_index+4] != '=O':
                                        if len(s) > co_index+2 and s[co_index+2] == '(': #  C(...)(...)OC(...
                                            index = s.find('(', co_index+2)
                                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                                            if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+2:co_index+4] != '=O':
                                                ether_number = ether_number + 1
                                                open_loop = False
                                                break

                                            else:
                                                i = i + 1

                                        else:
                                            ether_number = ether_number + 1
                                            open_loop = False
                                            break

                                    else:
                                        i = i + 1


                                else:
                                    i = i + 1
                    else:
                        # excluding ester:
                        if s[co_index+2:co_index+4] != '=O':
                            if len(s) > co_index+2 and s[co_index+2] == '(': #  C(...)(...)OC(...
                                index = s.find('(', co_index+2)
                                second_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O' and s[co_index+2:co_index+4] != '=O':
                                    ether_number = ether_number + 1


                            else:
                                ether_number = ether_number + 1


        co_index = s.find('OC', co_index + 1)


    # ether attached to a nonaromatic ring:
    if cyc_number != None:
        for j in range(1,cyc_number+1):
            A = s.find(str(j), 0)
            B = s.find(str(j), A+1)
            co_index = s.find('OC', 0)
            while co_index != -1:
                if s[co_index-2:co_index] == 'C'+str(j):  # C1OC
                    if B == co_index-1:  # 1...C1OC
                        # excluding ester:
                        if s[co_index+2:co_index+4] != '=O':
                            if len(s) > co_index+2 and s[co_index+2] == '(': #  ...C1OC(...
                                index = s.find('(', co_index+2)
                                first_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O':
                                    ether_number = ether_number + 1

                            else:
                                ether_number = ether_number + 1


                elif s[co_index-1] in [')']:  #  ...)OC...
                    index = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-2:first_opening_parenthesis] == 'C'+str(j):  #  ...C1(...)OC...
                        if B == first_opening_parenthesis-1 or first_opening_parenthesis < B < co_index-1:  #  ...C1(...1...)OC...
                            # excluding ester:
                            if s[co_index+2:co_index+4] != '=O':
                                if len(s) > co_index+2 and s[co_index+2] == '(': #  ...C1OC(...
                                    index = s.find('(', co_index+2)
                                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O':
                                        ether_number = ether_number + 1

                                else:
                                    ether_number = ether_number + 1


                elif s[co_index-1] in ['(']:  # ...(OC...
                    C = find_closing_parenthesis(s, co_index-1)
                    if s[co_index-3:co_index-1] == 'C'+str(j):  #  ...C1(OC...
                        if B == co_index-2 or B > C:  # 1...C1(OC...)...1
                            # excluding ester:
                            if s[co_index+2:co_index+4] != '=O':
                                if s[co_index+2] == '(': #  ...C1(OC(...
                                    index = s.find('(', co_index+2)
                                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O':
                                        ether_number = ether_number + 1

                                else:
                                    ether_number = ether_number + 1


                co_index = s.find('OC', co_index + 1)


    return ether_number

###################   Alicyclic ether  #######################################################
def alicyclic_ether(s):
    alicyclic_ether_number = 0
    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings

    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []


    O_pose = []   # list of 'O' index of alicyclic ethers in a string

    if cyc_number != None:
        for i in range(1,cyc_number+1):
            A = s.find(str(i), 0)
            while A != -1:  #  for the PRAM code in which rings have the same indexes
                B = s.find(str(i), A+1)
                if B !=-1:
                    C = s.find('OC', A, B+1)
                    while C != -1:
                        if s[C-1] not in ['O', 'N']:
                            D = s.find('(', A, B)
                            if D == -1:  #  ...C1...OC...C1...
                                if (C+2 == B) or ((C + 3 == B) and str(C + 2).isdigit()):  #  ...C1...COC1  or ...C2...COC12
                                    if s[C+3:C+5] == '=O':  #  ...C1...COC1=O
                                        C = -1
                                        pass

                                    elif (C-1 == A) or ((C-2 == A) and str(C-1).isdigit()):  #...C1OC1  or #...C12OC1
                                        if s[A-3:A-1] == 'O=':  # O=C1OC1
                                            C = -1
                                            pass

                                        else:  # C1OC1
                                            if C in O_pose:  #  check if O is not used for another cycle mistakenly
                                                C = -1
                                                pass

                                            else:
                                                O_pose.append(C)
                                                alicyclic_ether_number = alicyclic_ether_number + 1
                                                C = -1

                                    else:  #  ...C1...COC1
                                        if C in O_pose:
                                            C = -1
                                            pass

                                        else:
                                            O_pose.append(C)
                                            alicyclic_ether_number = alicyclic_ether_number + 1
                                            C = -1


                                elif (C-1 == A) or ((C-2 == A) and str(C-1).isdigit()):  #...C1OC....C1  or #...C12OC....C1
                                    if s[A-3:A-1] == 'O=':  # O=C1OC...C1
                                        C = s.find('OC', C+1, B+1)
                                        pass

                                    else:  # ...C1OC....C1
                                        if C in O_pose:
                                            C = s.find('OC', C+1, B+1)
                                            pass

                                        else:
                                            O_pose.append(C)
                                            alicyclic_ether_number = alicyclic_ether_number + 1
                                            C = s.find('OC', C+1, B+1)

                                else:  # C1...OC...C1
                                    if C in O_pose:
                                        C = s.find('OC', C+1, B+1)
                                        pass

                                    else:
                                        O_pose.append(C)
                                        alicyclic_ether_number = alicyclic_ether_number + 1
                                        C = s.find('OC', C+1, B+1)

                            else:
                                while D != -1:
                                    if D > C:  #  ...C1...OC..(...)...C1...
                                        if s[C+2:C+5] == '(=O':  #  ...C1...COC(=O)..(...)...C1...
                                            break

                                        elif (C-1 == A) or ((C-2 == A) and str(C-1).isdigit()):  #...C1OC....(...)...C1
                                            if s[A-3:A-1] == 'O=':  # O=C1OC...(...)...C1
                                                break

                                            else:  # ...C1OC....C1
                                                if C in O_pose:
                                                    break

                                                else:
                                                    O_pose.append(C)
                                                    alicyclic_ether_number = alicyclic_ether_number + 1
                                                    break

                                        else:
                                            if C in O_pose:
                                                break

                                            else:
                                                O_pose.append(C)
                                                alicyclic_ether_number = alicyclic_ether_number + 1
                                                break  # Break only from the inner while loop

                                    else:
                                        E = find_closing_parenthesis(s, D)
                                        if E > C:
                                            if E < B: #  ...C1...(...COC...)...C1...  ether is not alicyclic     [O]C1C(O)COC1=O
                                                break  # Break from the inner while loop

                                            elif s[C+2:C+5] == '(=O':  # ...C1...(...COC...C1...)...
                                                break

                                            else:
                                                if C in O_pose:
                                                    break

                                                else:
                                                    O_pose.append(C)
                                                    alicyclic_ether_number = alicyclic_ether_number + 1
                                                    break

                                        else:  #  ...C1...(...)...OC...C1...
                                            if s[C-4:C] == '(=O)':  # ...C1...(=O)OC...C1...
                                                break

                                            else:
                                                D = s.find('(', D + 1, B)
                                                if D == -1:  # No more parenthesis
                                                    if s[C+2:C+5] == '(=O':  # ...C1...(...)...COC(=O)...C1...
                                                        break

                                                    elif (C+2 == B) or ((C + 3 == B) and str(C + 2).isdigit()):  #  ...C1...COC1
                                                        if s[C+3:C+5] == '=O':  #  ...C1...COC1=O
                                                            break

                                                        else:
                                                            if C in O_pose:
                                                                break

                                                            else:
                                                                O_pose.append(C)
                                                                alicyclic_ether_number = alicyclic_ether_number + 1
                                                                break

                                                    else:
                                                        if C in O_pose:
                                                            break

                                                        else:
                                                            O_pose.append(C)
                                                            alicyclic_ether_number = alicyclic_ether_number + 1
                                                            break

                                C = s.find('OC', C+1, B+1)
                        else:
                            C = s.find('OC', C+1, B+1)

                else:
                    break



                A = s.find(str(i), B+1)


        for j in range(1,cyc_number+1):
            A = s.find(str(j), 0)
            while A != -1:
                if (s[A-1] == 'O') or (s[A-2] == 'O' and s[A-1].isdigit()):  #  O1 or O12
                    if len(s) > A+1 and s[A+1] in ['C', 'c'] :  #  O1C
                        B = s.find(str(j), A+1)
                        if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  O1C...C1 or O1C...C21
                            if s[A+2:A+6] != '(=O)' and s[B+1:B+3] != '=O':   #  O1C(=O)...C1  and  O1C...C1=O
                                alicyclic_ether_number = alicyclic_ether_number + 1

                    elif s[A-2] in ['C', 'c']:  #  CO1
                        B = s.find(str(j), 0, A)
                        if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C12...CO2
                            if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...CO1 or O=C1...CO1
                                alicyclic_ether_number = alicyclic_ether_number + 1

                    elif s[A-3:A-1] in list_tot:  #  C2O1
                        B = s.find(str(j), 0, A)
                        if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C2O1
                            if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...C2O1 or O=C1...C2O1
                                alicyclic_ether_number = alicyclic_ether_number + 1

                    elif s[A-2] == '(':  #  (O1
                        if s[A-3] in ['C', 'c']:  #  C(O1
                            B = s.find(str(j), 0, A)
                            if s[B-1] in ['C', 'c']:  #  C1...C(O1
                                if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  # O=C1...C(O1  or  C1(=O)...C(O1
                                    alicyclic_ether_number = alicyclic_ether_number + 1

                        elif s[A-4:A-2] in list_tot:  #  C2(O1
                            B = s.find(str(j), 0, A)
                            if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C2(O1
                                if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...C2(O1 or #  O=C1...C2(O1
                                    alicyclic_ether_number = alicyclic_ether_number + 1

                        elif s[A-3] == ')':  #  )(O1
                            index = s.rfind(')', 0, A-2)
                            first_opening_parenthesis = find_opening_parenthesis(s, index)
                            if s[first_opening_parenthesis-1] in ['C', 'c']:  #  C(...)(O1
                                B = s.find(str(j), 0, A)
                                if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C(...)(O1
                                    if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  O=C1...C(...)(O1  #  C1(=O)...C(...)(O1
                                        alicyclic_ether_number = alicyclic_ether_number + 1

                    elif s[A-2] == ')':  #  )O1
                        index = s.rfind(')', 0, A-1)
                        first_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[first_opening_parenthesis-1] in ['C', 'c']:  #  C(...)O1
                            B = s.find(str(j), 0, A)
                            if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C(...)O1
                                if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...C(...)O1  or  O=C1...C(...)O1
                                    if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  C1...C(=O)O1
                                        alicyclic_ether_number = alicyclic_ether_number + 1


                        elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C2(...)O1
                            B = s.find(str(j), 0, A)
                            if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C2(...)O1
                                if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...C2(...)O1  or  O=C1...C2(...)O1
                                    alicyclic_ether_number = alicyclic_ether_number + 1

                        elif s[first_opening_parenthesis-1] == ')':  #  )(...)O1
                            index = s.rfind(')', 0, first_opening_parenthesis)
                            second_opening_parenthesis = find_opening_parenthesis(s, index)
                            if s[second_opening_parenthesis-1] in ['C', 'c']:  #  C(...)(...)O1
                                B = s.find(str(j), 0, A)
                                if (s[B-1] in ['C', 'c']) or (s[B-2] in ['C', 'c'] and s[B-1].isdigit()):  #  C1...C(...)(...)O1
                                    if s[B+1:B+5] != '(=O)' or s[B-3:B-1] != 'O=':  #  C1(=O)...C(...)(...)O1  or  #  O=C1...C(...)(...)O1
                                        alicyclic_ether_number = alicyclic_ether_number + 1


                A = s.find(str(j), A+1)


    return alicyclic_ether_number

###################   Aromatic ether  #######################################################
def aromatic_ether_group(s):

    aromatic_ether_number = 0
    cyc_number = find_highest_digit(s)

    # Alkyl as the first part of the SMILE:
    co_index = s.find('Oc', 0)
    while co_index != -1:
        if s[co_index-1] in ['C'] and co_index-1 != -1:  #  COc...
            if s[co_index-3:co_index-1] != 'O=':  #  O=COc...
                aromatic_ether_number = aromatic_ether_number + 1

        elif s[co_index-1] in ['(']:  #  (Oc
            if s[co_index-2] in ['C']:  #  C(Oc
                index = s.find('(', co_index-1)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  # Avoid ester  C(Oc...)=O
                    aromatic_ether_number = aromatic_ether_number + 1

            elif s[co_index-2] in [')']:  #  )(Oc
                index = s.rfind(')', 0, co_index-1)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  #  C(...)(Oc
                    aromatic_ether_number = aromatic_ether_number + 1

        elif s[co_index-1] in [')']:  #  )Oc
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  C(...)Oc
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  C(=O)Oc...
                    aromatic_ether_number = aromatic_ether_number + 1

            elif s[first_opening_parenthesis-1] == ')':   #  )(...)Oc
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':   #  C(...)(...)Oc
                    aromatic_ether_number = aromatic_ether_number + 1

        co_index = s.find('Oc', co_index + 1)

    # ether attached to one aromatic and one non-aromatic rings:
    if cyc_number != None:
        for i in range(1,cyc_number+1):
            co_index = s.find('Oc', 0)
            while co_index != -1:
                if s[co_index-2:co_index] == 'C'+str(i) and co_index-2 != -2:  #  C1Oc
                    aromatic_ether_number = aromatic_ether_number + 1

                elif s[co_index-3:co_index] == 'C'+str(i) and co_index-3 != -3:  #  C12Oc
                    aromatic_ether_number = aromatic_ether_number + 1

                elif s[co_index-1] in ['(']:  #  (Oc
                    if s[co_index-3:co_index-1] == 'C'+str(i):  #  C1(Oc
                        aromatic_ether_number = aromatic_ether_number + 1

                    elif s[co_index-4:co_index-1] == 'C'+str(i):  #  C12(Oc
                        aromatic_ether_number = aromatic_ether_number + 1

                elif s[co_index-1] in [')']:  #  )Oc
                    index = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-2:first_opening_parenthesis] == 'C'+str(i):  #  C1(...)Oc
                        aromatic_ether_number = aromatic_ether_number + 1

                    elif s[first_opening_parenthesis-3:first_opening_parenthesis] == 'C'+str(i):  #  C12(...)Oc
                        aromatic_ether_number = aromatic_ether_number + 1

                co_index = s.find('Oc', co_index + 1)

    # Aromatic as the first part of the SMILE:
    if cyc_number != None:
        for i in range(1,cyc_number+1):
            co_index = s.find('OC', 0)
            while co_index != -1:
                if (s[co_index-2:co_index] == 'c'+str(i) and co_index-2 != -2) or (s[co_index-3:co_index] == 'c'+str(i) and co_index-3 != -3):  #  ...c1OC
                    if s[co_index+2:co_index+4] != '=O': #  ...c1OC=O/
                        if s[co_index+2] == '(': #  ...c1OC(...
                            index = s.find('(', co_index+2)
                            first_closing_parenthesis = find_closing_parenthesis(s, index)
                            if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...c1OC(=O)
                                aromatic_ether_number = aromatic_ether_number + 1

                        else:
                            aromatic_ether_number = aromatic_ether_number + 1

                elif s[co_index-1] in ['(']:  #  ...(OC
                    if s[co_index-3:co_index-1] == 'c'+str(i) or s[co_index-4:co_index-1] == 'c'+str(i):  #  ...c1(OC
                        if s[co_index+2:co_index+4] != '=O': #  ...c1(OC=O)...
                            if s[co_index+2] == '(': #  ...c1(OC(...
                                index = s.find('(', co_index+2)
                                first_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...c1(OC(=O)
                                    aromatic_ether_number = aromatic_ether_number + 1

                            else:
                                aromatic_ether_number = aromatic_ether_number + 1

                elif s[co_index-1] in [')']:  #  ...)OC
                    index = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-2:first_opening_parenthesis] == 'c'+str(i) or s[first_opening_parenthesis-3:first_opening_parenthesis] == 'c'+str(i):  #  c1(...)OC
                        if s[co_index+2:co_index+4] != '=O': #  ...c1(...)OC=O/
                            if s[co_index+2] == '(': #  ...c1(...)OC(...
                                index = s.find('(', co_index+2)
                                first_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...c1(...)OC(=O)
                                    aromatic_ether_number = aromatic_ether_number + 1

                            else:
                                aromatic_ether_number = aromatic_ether_number + 1

                co_index = s.find('OC', co_index + 1)

    co_index = s.find('OC', 0)
    while co_index != -1:
        if s[co_index-1] in ['(']:  #  ...(OC
            if s[co_index-2] in ['c']:  #  ...c(OC
                if s[co_index+2:co_index+4] != '=O': #  ...c(OC=O)
                    if s[co_index+2] == '(': #  ...c(OC(...
                        index = s.find('(', co_index+2)
                        first_closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...c(OC(=O)
                            aromatic_ether_number = aromatic_ether_number + 1

                    else:
                        aromatic_ether_number = aromatic_ether_number + 1

        elif s[co_index-1] in [')']:  #  ...)OC
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'c':  #  c(...)OC
                if s[co_index+2:co_index+4] != '=O': #  ...c(...)OC=O/
                    if len(s) > co_index+2 and s[co_index+2] == '(': #  ...c(...)OC(...
                        index = s.find('(', co_index+2)
                        first_closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O' and s[co_index+3:co_index+5] != '=O': #  ...c(...)OC(=O)
                            aromatic_ether_number = aromatic_ether_number + 1

                    else:
                        aromatic_ether_number = aromatic_ether_number + 1

        co_index = s.find('OC', co_index + 1)


    return aromatic_ether_number

###################   Nitrate  #######################################################
def nitrate_number(s):
    count_1 = s.count('ON(=O)=O')
    count_2 = s.count('O[N+](=O)[O-]')
    count_3 = s.count('O=N(=O)O')
    count_4 = s.count('O=[N+]([O-])O')
    count_5 = s.count('[O-][N+](=O)O')


    count_11 = s.count('OON(=O)=O')
    count_22 = s.count('OO[N+](=O)[O-]')
    count_33 = s.count('O=N(=O)OO')
    count_44 = s.count('O=[N+]([O-])OO')
    count_55 = s.count('[O-][N+](=O)OO')

    carbonylperoxynitrate = count_11 + count_22 + count_33 + count_44 + count_55


    nitrate = count_1 + count_2 + count_3 + count_4 + count_5

    return nitrate - carbonylperoxynitrate

###################   Nitro  #######################################################
def nitro_group(s):
    nitro_number = 0
    co_index = s.find('N(=O)=O', 0)
    while co_index != -1:
        if s[co_index-1] != 'O':
            nitro_number = nitro_number + 1
        co_index = s.find('N(=O)=O', co_index+1)

    co_index2 = s.find('O=N(=O)', 0)
    while co_index2 != -1:
        if s[co_index2+7] != 'O':
            nitro_number = nitro_number + 1
        co_index2 = s.find('O=N(=O)', co_index2+1)

    co_index = s.find('[N+](=O)[O-]', 0)
    while co_index != -1:
        if s[co_index-1] != 'O':
            nitro_number = nitro_number + 1
        co_index = s.find('[N+](=O)[O-]', co_index+1)

    co_index2 = s.find('O=[N+][O-]', 0)
    while co_index2 != -1:
        if s[co_index2+10] != 'O':
            nitro_number = nitro_number + 1
        co_index2 = s.find('O=[N+][O-]', co_index2+1)

    co_index3 = s.find('[O-][N+](=O)', 0)
    while co_index3 != -1:
        if s[co_index3+12] != 'O':
            nitro_number = nitro_number + 1
        co_index3 = s.find('[O-][N+](=O)', co_index3+1)


    return nitro_number

###################   Aromatic hydroxyl  #######################################################
def aromatic_hydroxyl_group(s):
    aromatic_hydroxyl_number = 0
    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_rings = []
    # hydroxyl as last characters of the SMILE
    if s[-1] == 'O':
        if s[-3:-1] in list_of_rings or s[-4:-1] in list_of_rings:  # ...c1O/ or # ...c12O/
            # excluding nitrophenol
            if s[-3:-1] in list_of_rings:
                index1 = s.rfind(s[-3:-1], 0, -3)
            else:
                index1 = s.rfind(s[-4:-1], 0, -3)

            if index1 != -1:

                if 'c(N(=O)=O)' in s[index1:-3]:  # c1...c(N(=O)=O)...c1O/
                    pass

                elif 'c([N+](=O)[O-])' in s[index1:-3]:  # c1...c(N(=O)=O)...c1O/
                    pass

                elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c1O/
                    pass

                elif s[-4:-1] in list_of_rings and s[index1+3:index1+12] == '(N(=O)=O)':  # c1(N(=O)=O)...c1O/
                    pass

                elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...c1O/
                    pass

                elif s[-4:-1] in list_of_rings and s[index1+3:index1+17] == '([N+](=O)[O-])':  # c1(N(=O)=O)...c1O/
                    pass

                elif s[index1-7:index1] == 'O=N(=O)':    # O=N(=O)c1...c1O/
                    pass

                elif s[index1-10:index1] == 'O=[N+][O-]':    # O=[N+][O-]c1...c1O/
                    pass

                elif s[index1-12:index1] == '[O-][N+](=O)':    # [O-][N+](=O)c1...c1O/
                    pass

                else:
                    aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1

        elif s[-2] == ')':  # ...)O
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            found = False   # Flag to control breaking both loops
            if s[first_opening_parenthesis-1] == 'c':  # ...c(...)O
                for j in list_of_rings:
                    if found:  # If the condition was met, break the for loop
                        break
                    index1 = s.find(j, 0, first_opening_parenthesis-1)
                    index2 = s.find(j, first_opening_parenthesis)
                    if index2 == -1:  # there is another cycle: c1...c1...c()O
                        pass

                    else:
                        index3 = s.find('c(', index1, index2)
                        while index3 != -1:
                            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                            if s[index3+2:index3+10] == 'N(=O)=O)': # c1...c(..c(N(=O)=O)...c1...)O
                                pass

                            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...c(..c([N+](=O)[O-])...c1...)O
                                pass

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...c(...c(...c1...)N(=O)=O)O
                                pass

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...c(...c(...c1...)[N+](=O)[O-])O
                                pass

                            elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c(..c1...)O
                                pass

                            elif len(j) == 3 and s[index1+3:index1+12] == '(N(=O)=O)':  # c1(N(=O)=O)...c(..c1...)O
                                pass

                            elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1([N+](=O)[O-])...c(..c1...)O
                                pass

                            elif len(j) == 3 and s[index1+3:index1+17] == '([N+](=O)[O-])':  # c1([N+](=O)[O-])...c(..c1...)O
                                pass

                            elif s[index2+2:index2+9] == 'N(=O)=O':  # c1...c(..c1N(=O)=O)O
                                pass

                            elif len(j) == 3 and s[index2+3:index2+10] == 'N(=O)=O':  # c1...c(..c1N(=O)=O)O
                                pass

                            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1...c(..c1[N+](=O)[O-])O  Cc1ccc(c(c1)C(C)(C)C)O
                                pass

                            elif len(j) == 3 and s[index2+3:index2+15] == '[N+](=O)[O-]':  # c1...c(..c1[N+](=O)[O-])O  Cc1ccc(c(c1)C(C)(C)C)O
                                pass

                            else:
                                if index3 == first_opening_parenthesis-1:  # avoid the hydroxyl conected carbon if there is another c( within its branch --> c(c(...)...)O
                                    index3 = s.find('c(', index3+1, index2)
                                    if index3 != -1:
                                        index3 = first_opening_parenthesis-1
                                        pass

                                    else:
                                        aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                        found = True
                                        break

                                else:
                                    aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                    found = True
                                    break

                            index3 = s.find('c(', index3+1, index2)

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings or s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  # c1(...)O
                if s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:
                    index1 = s[first_opening_parenthesis-2:first_opening_parenthesis]
                    index2 = s.find(s[first_opening_parenthesis-2:first_opening_parenthesis], first_opening_parenthesis)
                else:
                    index1 = s[first_opening_parenthesis-3:first_opening_parenthesis]
                    index2 = s.find(s[first_opening_parenthesis-3:first_opening_parenthesis], first_opening_parenthesis)

                if index2 != -1:
                    index3 = s.find('c(', index1, index2)
                    while index3 != -1:
                        first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                        if s[index3+2:index3+10] == 'N(=O)=O)': # c1(..c(N(=O)=O)...c1...)O
                            pass

                        elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1(..c(N(=O)=O)...c1...)O
                            pass

                        elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1(...c(...c1...)N(=O)=O)O
                            pass

                        elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1(...c(...c1...)N(=O)=O)O
                            pass

                        elif s[index2+2:index2+9] == 'N(=O)=O':  # c1(..c1N(=O)=O)O
                            pass

                        elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings and s[index2+3:index2+10] == 'N(=O)=O':  # c1(..c1N(=O)=O)O
                            pass

                        elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1(..c1N(=O)=O)O
                            pass

                        elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings and s[index2+3:index2+15] == '[N+](=O)[O-]':  # c1(..c1N(=O)=O)O
                            pass

                        else:
                            aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                            break

                        index3 = s.find('c(', index3+1, index2)


    #  hydroxyl as first characters of the SMILE
    if s.startswith('Oc1'):  # Oc1
        index2 = s.find('c1', 3)  # Oc1...c1
        index3 = s.find('c(', 2, index2)

        if index3 == -1:  # c1ccccc1
            aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1

        while index3 != -1:
            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
            if s[index3+2:index3+10] == 'N(=O)=O)': # Oc1...c(N(=O)=O)...c1...
                pass

            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # Oc1...c([N+](=O)[O-])...c1...
                pass

            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # Oc1...c(N(=O)=O...c1...)...
                pass

            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # Oc1...c([N+](=O)[O-]...c1...)...
                pass

            elif s[index2+2:index2+9] == 'N(=O)=O':  # Oc1...c1N(=O)=O
                pass

            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # Oc1...c1[N+](=O)[O-]
                pass

            else:
                aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                break

            index3 = s.find('c(', index3+1, index2)

    #  hydroxyl as middle characters of the SMILE
    co_index = s.find('c(O)', 0)
    while co_index != -1:  # ...c(O)...
        found = False
        for j in list_of_rings:
            if found:  # If the condition was met, break the for loop
                break
            index1 = s.find(j, 0, co_index)
            index2 = s.find(j, co_index)
            if index2 == -1:
                pass

            else:
                index3 = s.find('c(', index1, index2)
                while index3 != -1:
                    first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                    if s[index3+2:index3+10] == 'N(=O)=O)': # c1...c(N(=O)=O)...c(O)...c1...
                        pass

                    elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...c(N(=O)=O)...c(O)...c1...
                        pass

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...c(...c(O)...c1...)N(=O)=O...
                        pass

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...c(...c(O)...c1...)N(=O)=O...
                        pass

                    elif s[index2+2:index2+9] == 'N(=O)=O':  # c1...c(O)...c1N(=O)=O
                        pass

                    elif len(j) == 3 and s[index2+3:index2+10] == 'N(=O)=O':  # c1...c(O)...c12N(=O)=O
                        pass

                    elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1...c(O)...c1N(=O)=O
                        pass

                    elif len(j) == 3 and s[index2+3:index2+15] == '[N+](=O)[O-]':  # c12...c(O)...c12N(=O)=O
                        pass

                    elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c(O)...c1
                        pass

                    elif len(j) == 3 and s[index1+3:index1+12] == '(N(=O)=O)':  # c12(N(=O)=O)...c(O)...c1
                        pass

                    elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...c(O)...c1
                        pass

                    elif len(j) == 3 and s[index1+3:index1+17] == '([N+](=O)[O-])':  # c12(N(=O)=O)...c(O)...c1
                        pass

                    elif s[index1-7:index1] == 'O=N(=O)':  # O=N(=O)c1...c(O)...c1
                        pass

                    elif s[index1-10:index1] == 'O=[N+][O-]':  # O=N(=O)c1...c(O)...c1
                        pass

                    elif s[index1-12:index1] == '[O-][N+](=O)':  # [O-][N+](=O)c1...c(O)...c1
                        pass

                    else:
                        if index3 == co_index:  # avoid the hydroxyl conected carbon --> c(...)O
                            index3 = s.find('c(', index3+1, index2)
                            if index3 != -1:
                                index3 = co_index
                                pass

                            else:
                                aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                found = True
                                break

                        else:
                            aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                            found = True
                            break

                    index3 = s.find('c(', index3+1, index2)

        co_index = s.find('c(O)', co_index + 1)

    for j in list_of_rings:
        co_index2 = s.find(j + '(O)', 0)
        while co_index2 != -1:  # c1(O)...
            index = s.find(j, co_index2+4)
            index3 = s.find('c(', co_index2, index)
            while index3 != -1:
                first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                if s[index3+2:index3+10] == 'N(=O)=O)': # c1(O)...c(N(=O)=O)...c1...
                    pass

                elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1(O)...c([N+](=O)[O-])...c1...
                    pass

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1(O)...c(...c1...)N(=O)=O...
                    pass

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1(O)...c(...c1...)[N+](=O)[O-]...
                    pass

                elif s[index2+2:index2+9] == 'N(=O)=O':  # c1(O)...c1N(=O)=O
                    pass

                elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1(O)...c1N(=O)=O
                    pass

                else:
                    aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                    break

                index3 = s.find('c(', index3+1, index)


            co_index2 = s.find(j + '(O)', co_index2 + 1)

    #  hydroxyl as last characters of a branch in the SMILE
    co_index3 = s.find('O)', 0)
    while co_index3 != -1:
        if s[co_index3-2:co_index3] in list_of_rings or s[co_index3-3:co_index3] in list_of_rings:  # ...(...c12O)
            if  s[co_index3-2:co_index3] in list_of_rings:
                index1 = s.rfind(s[co_index3-2:co_index3], 0, co_index3-2)

            else:
                index1 = s.rfind(s[co_index3-3:co_index3], 0, co_index3-2)

            if index1 != -1:

                index3 = s.find('c(', index1, co_index3-2)

                while index3 != -1:
                    first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                    if s[index3+2:index3+10] == 'N(=O)=O)': # c1...(...c(N(=O)=O)...c1O)
                        pass

                    elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...(...c(N(=O)=O)...c1O)
                        pass

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...(...c(...c1O...)N(=O)=O)
                        pass

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...(...c(...c1O...)N(=O)=O)
                        pass

                    elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...(...c1O)
                        pass

                    elif s[co_index3-3:co_index3] in list_of_rings and s[index1+3:index1+12] == '(N(=O)=O)':  # c1(N(=O)=O)...(...c1O)
                        pass

                    elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...(...c1O)
                        pass

                    elif s[co_index3-3:co_index3] in list_of_rings and s[index1+3:index1+17] == '([N+](=O)[O-])':  # c1(N(=O)=O)...(...c1O)
                        pass

                    else:
                        aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                        break

                    index3 = s.find('c(', index3+1, co_index3-2)


        elif s[co_index3-1] == ')':  # ...)O)
            index = s.rfind(')', 0, co_index3)
            second_opening_parenthesis = find_opening_parenthesis(s, index)
            index_0 = s.rfind(')', 0, co_index3+2)
            first_opening_parenthesis = find_opening_parenthesis(s, index_0)

            if s[second_opening_parenthesis-1] == 'c':  # ...c(...)O)...
                found = False
                for j in list_of_rings:
                    if found:  # If the condition was met, break the for loop
                        break
                    index1 = s.find(j, 0, second_opening_parenthesis-1)
                    index2 = s.find(j, second_opening_parenthesis)
                    if index2 == -1:
                        pass

                    else:
                        index3 = s.find('c(', index1, index2)
                        while index3 != -1:
                            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                            if s[index3+2:index3+10] == 'N(=O)=O)': # (...c1...c(N(=O)=O)...c(...c1...)O)...
                                pass

                            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # (...c1...c([N+](=O)[O-])...c(...c1...)O)...
                                pass

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # (...c1...c(...c(...c1...)N(=O)=O)O)...
                                pass

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # (...c1...c(...c(...c1...)N(=O)=O)O)...
                                pass

                            elif s[index2+2:index2+9] == 'N(=O)=O':  # (...c1...c(..c1N(=O)=O)O
                                pass

                            elif len(j) == 3 and s[index2+3:index2+10] == 'N(=O)=O':  # (...c1...c(..c1N(=O)=O)O
                                pass

                            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # (...c1...c(..c1N(=O)=O)O
                                pass

                            elif len(j) == 3 and s[index2+3:index2+15] == '[N+](=O)[O-]':  # (...c1...c(..c1N(=O)=O)O
                                pass

                            else:   # c1cc(cc(c1)O)[N+](=O)[O-]
                                if index3 == second_opening_parenthesis-1:  # avoid the hydroxyl conected carbon --> c(...)O
                                    index3 = s.find('c(', index3+1, index2)
                                    if index3 != -1:
                                        index3 = second_opening_parenthesis-1
                                        pass

                                    else:
                                        if s[first_opening_parenthesis-1] == 'c' or s[first_opening_parenthesis-2:first_opening_parenthesis] == j:
                                            if s[co_index3+2:co_index3+9] == 'N(=O)=O': # ...c1...c(...c(c1...)O)N(=O)=O
                                                found = True
                                                break

                                            elif s[co_index3+2:co_index3+14] == '[N+](=O)[O-]': # ...c1...c(...c(c1...)O)[N+](=O)[O-]
                                                found = True
                                                break

                                            else:
                                                aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                                found = True
                                                break

                                        else:
                                            aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                            found = True
                                            break


                                else:
                                    aromatic_hydroxyl_number = aromatic_hydroxyl_number + 1
                                    found = True
                                    break

                            index3 = s.find('c(', index3+1, index2)


        co_index3 = s.find('O)', co_index3 + 1)


    return aromatic_hydroxyl_number

###################   Primary_amine  #######################################################
def primary_amine_group(s):
    primary_amine_number = 0
    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_rings = ['C' + str(i) for i in range(1, cyc_number+1)]

    else:
        list_of_rings = []

    # primary_amine as last characters of the SMILE:
    if s[-1] == 'N':  #  N/
        if s[-2] == 'C':  #  CN/
            if s[-4:-2] != 'O=':  #  O=CN
                primary_amine_number = primary_amine_number + 1

        elif s[-2] == ')':  #  )N/
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  C(...)N/
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  C(=O)N/
                    if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=' and first_opening_parenthesis-1 != 0:  #  O=C(...)N/ and /C(...)N/
                        primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  )(...)N/
                index2 = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index2)
                if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)N/
                    if s[second_opening_parenthesis+1:second_opening_parenthesis+3] != '=O':  #  C(=O)(...)N/
                        if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O' and second_opening_parenthesis-1 != 0:  #  C(...)(=O)N/ and /C(...)(...)N/
                            primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings and first_opening_parenthesis-2 != 0:  #  C1(...)N/  and  /C1(...)N/
                primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  #  C12(...)N/  and  /C1(...)N/
                primary_amine_number = primary_amine_number + 1

        elif s[-3:-1] in list_of_rings:  #  C1N/
            primary_amine_number = primary_amine_number + 1

        elif s[-4:-1] in list_of_rings:  #  C12N/
            primary_amine_number = primary_amine_number + 1


    #  primary_amine as first characters of the SMILE:
    #  excluding if the string starts with primary_amine i.e., /CN/ since they were included above
    if s.startswith('C('):  #  /C(
        if s[2:4] != '=O':  #  /C(=O)
            if s[2:4] == 'N)':  #  /C(N)
                if s[4:6] != '=O':  #  /C(N)=O
                    if s[4] == '(':  #  /C(N)(
                        if s[5:7] != '=O':  #  /C(N)(=O)
                            if s[5:7] == 'N)':  #  /C(N)(N)
                                if s[7:9] != '=O':  #  /C(N)(N)=O
                                    if s[7] == 'N' and len(s) == 8:  #  /C(N)(N)N/
                                        primary_amine_number = primary_amine_number + 3

                                    elif s[7] == '(':  #  /C(N)(N)(
                                        if s[8:10] == 'N)':  #  /C(N)(N)(N)
                                            if s[10] == 'N' and len(s) == 11:  #  /C(N)(N)(N)N/
                                                primary_amine_number = primary_amine_number + 4

                                            else:  #  /C(N)(N)(N)...
                                                primary_amine_number = primary_amine_number + 3

                                        else:
                                            index = s.find('(', 7)
                                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                                            if s[second_closing_parenthesis+1] == 'N' and len(s) == second_closing_parenthesis+2:  #  /C(N)(N)(...)N/
                                                primary_amine_number = primary_amine_number + 3

                                            else:  #  /C(N)(N)(...)...
                                                primary_amine_number = primary_amine_number + 2

                                    else:  #  /C(N)(N)...
                                        primary_amine_number = primary_amine_number + 2

                            else:
                                index = s.find('(', 4)
                                first_closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  /C(N)(...)=O
                                    if s[first_closing_parenthesis+1] == 'N' and len(s) == first_closing_parenthesis+2:  #  /C(N)(...)N/
                                        primary_amine_number = primary_amine_number + 2

                                    elif s[first_closing_parenthesis+1] == '(':  #  /C(N)(...)(
                                        if s[first_closing_parenthesis+2:first_closing_parenthesis+4] == 'N)':  #  /C(N)(...)(N)
                                            if s[first_closing_parenthesis+4] == 'N' and len(s) == first_closing_parenthesis+5:  #  /C(N)(...)(N)N/
                                                primary_amine_number = primary_amine_number + 3

                                            else:  #  /C(N)(...)(N)...
                                                primary_amine_number = primary_amine_number + 2


                                        else:
                                            index = s.find('(', first_closing_parenthesis+1)
                                            second_closing_parenthesis = find_closing_parenthesis(s, index)
                                            if s[second_closing_parenthesis+1] == 'N' and len(s) == second_closing_parenthesis+2:  #  /C(N)(...)(...)N/
                                                primary_amine_number = primary_amine_number + 2

                                            else:  #  /C(N)(...)(...)...
                                                primary_amine_number = primary_amine_number + 1

                                    else:  #  /C(N)(...)...
                                        primary_amine_number = primary_amine_number + 1

                    elif s[4] == 'N' and len(s) == 5:  #  /C(N)N/
                        primary_amine_number = primary_amine_number + 2

                    else:  #  /C(N)...
                        primary_amine_number = primary_amine_number + 1

            else:
                index = s.find('(', 1)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] != '=O':  #  /C(...)=O
                    if s[first_closing_parenthesis+1] == 'N' and len(s) == first_closing_parenthesis+2:  #  /C(...)N/
                        primary_amine_number = primary_amine_number + 1

                    elif s[first_closing_parenthesis+1] == '(':  #  /C(...)(
                        if s[first_closing_parenthesis+2:first_closing_parenthesis+4] == 'N)':  #  /C(...)(N)...
                            if s[first_closing_parenthesis+4:first_closing_parenthesis+6] != '=O':  #  /C(...)(N)=O
                                if s[first_closing_parenthesis+4] == 'N' and len(s) == first_closing_parenthesis+5:  #  /C(...)(N)N/
                                    primary_amine_number = primary_amine_number + 2

                                elif s[first_closing_parenthesis+4] == '(':  #  /C(...)(N)(
                                    if s[first_closing_parenthesis+5:first_closing_parenthesis+7] == 'N)':  #  /C(...)(N)(N)
                                        if s[first_closing_parenthesis+7] == 'N' and len(s) == first_closing_parenthesis+8:  #  /C(...)(N)(N)N/
                                            primary_amine_number = primary_amine_number + 3

                                        else: #  /C(...)(N)(N)...
                                            primary_amine_number = primary_amine_number + 2

                                    else:
                                        index = s.find('(', first_closing_parenthesis+4)
                                        second_closing_parenthesis = find_closing_parenthesis(s, index)
                                        if s[second_closing_parenthesis+1] == 'N' and len(s) == second_closing_parenthesis+2:  #  /C(...)(N)(...)N/
                                            primary_amine_number = primary_amine_number + 2

                                        else:  #  /C(...)(N)(...)...
                                            primary_amine_number = primary_amine_number + 1


                                else:  #  /C(...)(N)...
                                    primary_amine_number = primary_amine_number + 1


    if s.startswith('NC'):  #  /NC
        if len(s) > 2 and s[2] != '1':  #  /NC1
            if s[2:4] != '=O':  #  /NC=O
                if s[2] == 'N':  #  /NCN
                    primary_amine_number = primary_amine_number + 2  # /NC(N)... can have two primary_amine groups considered latter

                else:  #  /NC...
                    primary_amine_number = primary_amine_number + 1

    # nonaromatic ring-connected primary_amine
    for j in list_of_rings:
        if s.startswith(j + '('):  #  /C1(
            if s[3:5] == 'N)':  #  /C1(N)
                if s[5] == '(':  #  /C1(N)(
                    if s[6:8] == 'N)':  #  /C1(N)(N)
                        primary_amine_number = primary_amine_number + 2

                    else:
                        primary_amine_number = primary_amine_number + 1

                elif s[5:7] != '=O':  #  /C1(N)=O
                    primary_amine_number = primary_amine_number + 1

            elif s[3:5] != '=O':  #  /C1(=O)
                index = s.find('(', 2)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '(N)':  #  /C1(...)(N)...
                    if s[first_closing_parenthesis+4] == 'N':  #  /C1(...)(N)N
                        primary_amine_number = primary_amine_number + 2

                    else:
                        primary_amine_number = primary_amine_number + 1

        if s.startswith('N' + j):  #  /NC1
            if s[3] == '(':  #  /NC1(
                if s[4:6] == 'N)':   #  /NC1(N)
                    primary_amine_number = primary_amine_number + 2

                else:
                    primary_amine_number = primary_amine_number + 1

            else:
                primary_amine_number = primary_amine_number + 1


    #  primary_amine as middle characters of the SMILE:
    co_index = s.find('(N)', 0)
    while co_index != -1:
        if s[co_index-1] == 'C':  #  C(N)
            if s[co_index-3:co_index-1] != 'O=':  #  O=C(N)
                if s[co_index+3:co_index+5] != '=O':  #  C(N)=O
                    if co_index-1 != 0:  #  /C(N)!
                        primary_amine_number = primary_amine_number + 1

        elif s[co_index-2:co_index] in list_of_rings:  #  C1(N)...
            if co_index-2 != 0:  #  /C1(N)...!
                primary_amine_number = primary_amine_number + 1

        elif s[co_index-3:co_index] in list_of_rings:  #  C12(N)...
            primary_amine_number = primary_amine_number + 1

        elif s[co_index-1] == ')':  #  )(N)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  C(...)(N)...
                if first_opening_parenthesis-1 != 0:  #  /C(...)(N)...!
                    primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  #  C1(...)(N)...
                if first_opening_parenthesis-2 != 0:  #  /C1(...)(N)...!
                    primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  #  C12(...)(N)...
                primary_amine_number = primary_amine_number + 1

        co_index = s.find('(N)', co_index + 1)


    #  primary_amine as last characters of a branch in the SMILE:
    co_index2 = s.find('N)', 0)
    while co_index2 != -1:
        if s[co_index2-1] == 'C':  # CN)
            primary_amine_number = primary_amine_number + 1

        elif s[co_index2-2:co_index2] in list_of_rings:  # C1N)
            primary_amine_number = primary_amine_number + 1

        elif s[co_index2-3:co_index2] in list_of_rings:  # C12N)
            primary_amine_number = primary_amine_number + 1

        elif s[co_index2-1] == ')':  # )N)
            index = s.rfind(')', 0, co_index2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  # ...C(...)N)
                primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  # ...C1(...)N)
                primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  # ...C12(...)N)
                primary_amine_number = primary_amine_number + 1

            elif s[first_opening_parenthesis-1] == ')':  # ...)(...)N)
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  # C(...)(...)N)
                    primary_amine_number = primary_amine_number + 1

        co_index2 = s.find('N)', co_index2 + 1)


    return primary_amine_number

###################   Secondary_amine  #######################################################
def secondary_amine_group(s):
    secondary_amine_number = 0
    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_of_nonarom_rings = list_of_nonarom_rings + list_of_arom_rings
    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_of_nonarom_rings = []
    # secondary amine as last characters of the SMILE
    if s[-2:] == 'NC':  #  NC/
        if len(s) >= 3 and s[-3] in ['C']:  #  CNC/
            if len(s) == 3:  #  /CNC/
                secondary_amine_number = secondary_amine_number + 1

            elif s[-5:-3] != 'O=':
                secondary_amine_number = secondary_amine_number + 1

        elif len(s) >= 3 and s[-3] == ')':  #  )NC
            index = s.rfind(')', 0, -2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C']:  #  C(...)NC
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':
                    if first_opening_parenthesis-1 == 0:  #  /C(...)NC
                        secondary_amine_number = secondary_amine_number + 1

                    elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':
                        secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_nonarom_rings:  #  C1(...)NC
                secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_nonarom_rings:  #  C12(...)NC
                secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  )(...)NC
                index1 = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index1)
                if s[second_opening_parenthesis-1] in ['C']:  #  C(...)(...)NC
                    secondary_amine_number = secondary_amine_number + 1

        elif s[-4:-2] in list_of_nonarom_rings:  #  C1NC/
            secondary_amine_number = secondary_amine_number + 1

        elif s[-5:-2] in list_of_nonarom_rings:  #  C12NC/
            secondary_amine_number = secondary_amine_number + 1

    # secondary amine as middle characters of the SMILE:
    co_index = s.find('N', 0)

    if len(s) >= co_index+2 and s[co_index+1] in ['C']:  #  NC or Nc
        while co_index != -1:
            if s[co_index-1] in ['C'] and co_index-1 != -1:  #  ...CNC
                if s[co_index-3:co_index-1] != 'O=':
                    if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                        if len(s) > co_index+2 and s[co_index+2] == '(':  #  ...CNC(
                            index = s.find('(', co_index+2)
                            closing_parenthesis = find_closing_parenthesis(s, index)
                            if s[co_index+3:co_index+5] != '=O':
                                if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                    secondary_amine_number = secondary_amine_number + 1

                        else:
                            secondary_amine_number = secondary_amine_number + 1

            elif s[co_index-2:co_index] in list_of_nonarom_rings:  #  ...C1NC
                if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                    if len(s) > co_index+2 and s[co_index+2] == '(':  #  ...C1NC(
                        index = s.find('(', co_index+2)
                        closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[co_index+3:co_index+5] != '=O':
                            if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                secondary_amine_number = secondary_amine_number + 1

                    else:
                        secondary_amine_number = secondary_amine_number + 1


            elif s[co_index-3:co_index] in list_of_nonarom_rings:  #  ...C12NC
                if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                    if len(s) > co_index+2 and s[co_index+2] == '(':  #  ...C12NC(
                        index = s.find('(', co_index+2)
                        closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[co_index+3:co_index+5] != '=O':
                            if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                secondary_amine_number = secondary_amine_number + 1

                    else:
                        secondary_amine_number = secondary_amine_number + 1

            elif s[co_index-1] == ')':  #  ...)NC
                index = s.rfind(')', 0, co_index)
                opening_parenthesis = find_opening_parenthesis(s, index)
                if s[opening_parenthesis-1] in ['C']:  #  C(...)NC
                    if s[opening_parenthesis+1:opening_parenthesis+3] != 'O=':
                        if s[opening_parenthesis-3:opening_parenthesis-1] != 'O=':
                            if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                                if len(s) > co_index+2 and s[co_index+2] == '(':  #  C(...)NC(...
                                    index = s.find('(', co_index+2)
                                    closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[co_index+3:co_index+5] != '=O':
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                            secondary_amine_number = secondary_amine_number + 1

                                else:
                                    secondary_amine_number = secondary_amine_number + 1

                elif s[opening_parenthesis-2:opening_parenthesis] in list_of_nonarom_rings:  #  C1(...)NC
                    if s[opening_parenthesis+1:opening_parenthesis+3] != 'O=':  # No need but anyway...
                        if s[opening_parenthesis-3:opening_parenthesis-1] != 'O=':
                            if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                                if len(s) > co_index+2 and s[co_index+2] == '(':  #  C1(...)NC(...
                                    index = s.find('(', co_index+2)
                                    closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[co_index+3:co_index+5] != '=O':
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                            secondary_amine_number = secondary_amine_number + 1

                                else:
                                    secondary_amine_number = secondary_amine_number + 1


                elif s[opening_parenthesis-3:opening_parenthesis] in list_of_nonarom_rings:  #  C12(...)NC
                    if s[opening_parenthesis+1:opening_parenthesis+3] != 'O=':  # No need but anyway...
                        if s[opening_parenthesis-3:opening_parenthesis-1] != 'O=':
                            if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                                if len(s) > co_index+2 and s[co_index+2] == '(':  #  C12(...)NC(...
                                    index = s.find('(', co_index+2)
                                    closing_parenthesis = find_closing_parenthesis(s, index)
                                    if s[co_index+3:co_index+5] != '=O':
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                            secondary_amine_number = secondary_amine_number + 1

                                else:
                                    secondary_amine_number = secondary_amine_number + 1

                elif s[opening_parenthesis-1] == ')':  #  ...)(...)NC
                    index = s.rfind(')', 0, opening_parenthesis)
                    opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[opening_parenthesis-1] in ['C']:  #  C(...)(...)NC
                        if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                            if len(s) > co_index+2 and s[co_index+2] == '(':  #  C(...)(...)NC(
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                    if s[closing_parenthesis+1] == '(':
                                        secondary_amine_number = secondary_amine_number + 1

                                    else:
                                        if s[co_index+3:co_index+5] != '=O':
                                            secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1

            elif s[co_index-1] == '(' and co_index-1 != -1:  #  ...(NC
                if s[co_index-2] in ['C']:  #  C(NC
                    index = s.find('(', co_index-1)
                    closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[co_index-4:co_index-2] != 'O=' and s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                        if s[co_index+2:co_index+4] != '=O':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1

                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1


                elif s[co_index-3:co_index-1] in list_of_nonarom_rings:  #  C1(NC
                    index = s.find('(', co_index-1)
                    closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[co_index-4:co_index-2] != 'O=' and s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                        if s[co_index+2:co_index+4] != '=O':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1

                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1

                elif s[co_index-4:co_index-1] in list_of_nonarom_rings:  #  C12(NC
                    index = s.find('(', co_index-1)
                    closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[co_index-4:co_index-2] != 'O=' and s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                        if s[co_index+2:co_index+4] != '=O':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1

                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1

                elif s[co_index-2] == ')':  #  ...)(NC
                    index = s.rfind(')', 0, co_index-1)
                    opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[opening_parenthesis-1] in ['C']:  #  C(...)(NC
                        if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1
                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1

                    elif s[opening_parenthesis-2:opening_parenthesis] in list_of_nonarom_rings:  #  C1(...)(NC
                        if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1
                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1

                    elif s[opening_parenthesis-3:opening_parenthesis] in list_of_nonarom_rings:  #  C12(...)(NC
                        if s[co_index+2:co_index+4] != '=O' and co_index+1 != len(s)-1 and s[co_index+2] != ')':
                            if s[co_index+2] == ')':
                                secondary_amine_number = secondary_amine_number + 1
                            elif s[co_index+2] == '(':
                                index = s.find('(', co_index+2)
                                closing_parenthesis = find_closing_parenthesis(s, index)
                                if s[co_index+3:co_index+5] != '=O':
                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':
                                        secondary_amine_number = secondary_amine_number + 1

                            else:
                                secondary_amine_number = secondary_amine_number + 1


            co_index = s.find('NC', co_index + 1)

    # secondary amine as last characters of a branch in the SMILE
    co_index = s.find('NC)', 0)
    while co_index != -1:
        if s[co_index-1] in ['C']:  #  ...CNC)
            if co_index-1 == 0:
                secondary_amine_number = secondary_amine_number + 1

            elif s[co_index-3:co_index-1] != 'O=':
                secondary_amine_number = secondary_amine_number + 1

        elif s[co_index-2:co_index] in list_of_nonarom_rings:  #  ...C1NC)
            secondary_amine_number = secondary_amine_number + 1

        elif s[co_index-3:co_index] in list_of_nonarom_rings:  #  ...C12NC)
            secondary_amine_number = secondary_amine_number + 1

        elif s[co_index-1] == ')':  #  ...)NC)
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C']:  #  C(...)NC)
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':
                    if first_opening_parenthesis-1 == 0:
                        secondary_amine_number = secondary_amine_number + 1

                    elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':
                        secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_nonarom_rings:  #  C1(...)NC)
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':
                    secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_nonarom_rings:  #  C12(...)NC)
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':
                    secondary_amine_number = secondary_amine_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  ...)(...)NC)
                index1 = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index1)
                if s[second_opening_parenthesis-1] in ['C']:  #  C(...)(...)NC)
                    secondary_amine_number = secondary_amine_number + 1

        co_index = s.find('NC)', co_index + 1)


    co_index = s.find('N', 0)
    while co_index != -1:
        if s[co_index-1] == 'C' and co_index-1 != -1:  #  ...CN
            if len(s) >= co_index+2 and s[co_index+1] == '(':  #  ...CN(
                if s[co_index+2] == 'C':  #  ...CN(C
                    if s[co_index+3:co_index+5] != '=O':  #  ...CN(C=O)
                        if s[co_index+3] == ')':  #  ...CN(C)
                            if s[co_index+4] == 'O':  #  ...CN(C)O
                                if s[co_index-3:co_index-1] != 'O=':  #  O=CN(C)O
                                    secondary_amine_number = secondary_amine_number + 1


                        elif s[co_index+3] == '(':  #  ...CN(C(
                            if s[co_index+4:co_index+6] != '=O':  #  ...CN(C(=O)
                                closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...CN(C(A)=O
                                    second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[second_closing_parenthesis+1] == 'O':  #  ...CN(C(A)B)O
                                        if s[co_index-3:co_index-1] != 'O=':  #  O=CN(C(A)B)O
                                            secondary_amine_number = secondary_amine_number + 1


                        else:
                            closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                            if s[closing_parenthesis+1] == 'O':  #  ...CN(CA)O
                                if s[co_index-3:co_index-1] != 'O=':  #  O=CN(CA)O
                                    secondary_amine_number = secondary_amine_number + 1


        elif s[co_index-1] == ')' and co_index-1 != -1:  #  ...)N...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  ...C(R)N...
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  ...C(=O)N...

                    if len(s) >= co_index+2 and s[co_index+1] == '(':  #  ...C(R)N(
                        if s[co_index+2] == 'C':  #  ...C(R)N(C
                            if s[co_index+3:co_index+5] != '=O':  #  ...C(R)N(C=O)
                                if s[co_index+3] == ')':  #  ...C(R)N(C)
                                    if s[co_index+4] == 'O':  #  ...C(R)N(C)O
                                        if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(C)C
                                            secondary_amine_number = secondary_amine_number + 1


                                elif s[co_index+3] == '(':  #  ...C(R)N(C(
                                    if s[co_index+4:co_index+6] != '=O':  #  ...C(R)N(C(=O)
                                        closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)N(C(A)=O
                                            second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                            if s[second_closing_parenthesis+1] == 'O':  #  ...C(R)N(C(A)B)O
                                                if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(C(A)B)O
                                                    secondary_amine_number = secondary_amine_number + 1


                                else:
                                    closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[closing_parenthesis+1] == 'O':  #  ...C(R)N(CA)O
                                        if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(CA)O
                                            secondary_amine_number = secondary_amine_number + 1


            elif s[first_opening_parenthesis-1] == ')':  #  ...)(R)N...
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  #  ...C(R)(R')N...

                    if len(s) >= co_index+2 and s[co_index+1] == '(':  #  ...C(R)(R')N(
                        if s[co_index+2] == 'C':  #  ...C(R)(R')N(C
                            if s[co_index+3:co_index+5] != '=O':  #  ...C(R)(R')N(C=O)
                                if s[co_index+3] == ')':  #  ...C(R)(R')N(C)
                                    if s[co_index+4] == 'O':  #  ...C(R)(R')N(C)O
                                        secondary_amine_number = secondary_amine_number + 1


                                elif s[co_index+3] == '(':  #  ...C(R)(R')N(C(
                                    if s[co_index+4:co_index+6] != '=O':  #  ...C(R)(R')N(C(=O)
                                        closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)(R')N(C(A)=O
                                            second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                            if s[second_closing_parenthesis+1] == 'O':  #  ...C(R)(R')N(C(A)B)O
                                                secondary_amine_number = secondary_amine_number + 1


                                else:
                                    closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[closing_parenthesis+1] == 'O':  #  ...C(R)(R')N(CA)O
                                        secondary_amine_number = secondary_amine_number + 1


        co_index = s.find('N', co_index + 1)


    if cyc_number != None: # if there is any cycle
    # Ending with as a cycle
        for i in range(1, cyc_number+1):
            co_index = s.find('(N' + str(i) + ')', 0)  #   ...(N1)
            while co_index != -1:
                if s[co_index-1] == 'C':  # ...C(N1)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...C(N1)...
                        secondary_amine_number = secondary_amine_number + 1

                elif s[co_index-1] == ')':  # ...)(N1)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)(N1)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)(N1)...
                            secondary_amine_number = secondary_amine_number + 1

                co_index = s.find('(N' + str(i) + ')', co_index + 1)

        for i in range(1, cyc_number+1):
            end_str = 'N' + str(i)
            if s.endswith(end_str):
                if len(s) > len(end_str) and s[-len(end_str) - 1] == 'C':   # CN1/
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...CN1/
                        secondary_amine_number = secondary_amine_number + 1

                elif len(s) > len(end_str) and s[-len(end_str) - 1] == ')':   # )N1/
                    index1 = s.rfind(')', 0, -len(end_str))
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)N1/
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)N1/
                            secondary_amine_number = secondary_amine_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)N1/
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)N1/
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)N1/
                                secondary_amine_number = secondary_amine_number + 1


    return secondary_amine_number

###################   Tertiary_amine  #######################################################
def tertiary_amine_group(s):
    tertiary_amine_number = 0
    co_index = s.find('N', 0)
    while co_index != -1:
        if s[co_index-1] == 'C' and co_index-1 != -1:  #  ...CN
            if len(s) > co_index+1 and s[co_index+1] == '(':  #  ...CN(
                if s[co_index+2] == 'C':  #  ...CN(C
                    if s[co_index+3:co_index+5] != '=O':  #  ...CN(C=O)
                        if s[co_index+3] == ')':  #  ...CN(C)
                            if s[co_index+4] == 'C':  #  ...CN(C)C
                                if s[co_index-3:co_index-1] != 'O=':  #  O=CN(C)C
                                    if s[co_index+5:co_index+7] != '=O':  #  CN(C)C=O
                                        if len(s) > co_index+5 and s[co_index+5] == '(':  #  CN(C)C(
                                            if s[co_index+6:co_index+8] != '=O':  #  CN(C)C(=O)
                                                closing_parenthesis = find_closing_parenthesis(s, co_index+5)
                                                if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  CN(C)C(...)=O
                                                    tertiary_amine_number = tertiary_amine_number + 1

                                        else:
                                            tertiary_amine_number = tertiary_amine_number + 1


                        elif s[co_index+3] == '(':  #  ...CN(C(
                            if s[co_index+4:co_index+6] != '=O':  #  ...CN(C(=O)
                                closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...CN(C(A)=O
                                    second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[second_closing_parenthesis+1] == 'C':  #  ...CN(C(A)B)C
                                        if s[co_index-3:co_index-1] != 'O=':  #  O=CN(C(A)B)C
                                            if s[second_closing_parenthesis+2:second_closing_parenthesis+4] != '=O':  #  ...CN(C(A)B)C=O
                                                if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+2] == '(':  #  ...CN(C(A)B)C(
                                                    if s[second_closing_parenthesis+3:second_closing_parenthesis+5] != '=O':  #  ...CN(C(A)B)C(=O)
                                                        third_closing_parenthesis = find_closing_parenthesis(s, second_closing_parenthesis+2)
                                                        if s[third_closing_parenthesis+1:third_closing_parenthesis+3] != '=O':  #  ...CN(C(A)B)C(I)=O
                                                            tertiary_amine_number = tertiary_amine_number + 1

                                                else:
                                                    tertiary_amine_number = tertiary_amine_number + 1


                        else:
                            closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                            if s[closing_parenthesis+1] == 'C':  #  ...CN(CA)C
                                if s[co_index-3:co_index-1] != 'O=':  #  O=CN(CA)C
                                    if s[closing_parenthesis+2:closing_parenthesis+4] != '=O':  #  CN(CA)C=O
                                        if len(s) > closing_parenthesis+2 and s[closing_parenthesis+2] == '(':  #  CN(CA)C(
                                            if s[closing_parenthesis+3:closing_parenthesis+5] != '=O':  #  CN(CA)C(=O)
                                                second_closing_parenthesis = find_closing_parenthesis(s, closing_parenthesis+2)
                                                if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O':  #  CN(CA)C(...)=O
                                                    tertiary_amine_number = tertiary_amine_number + 1

                                        else:
                                            tertiary_amine_number = tertiary_amine_number + 1



        elif s[co_index-1] == ')':  #  ...)N...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  ...C(R)N...
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  ...C(=O)N...

                    if len(s) > co_index+1 and s[co_index+1] == '(':  #  ...C(R)N(
                        if s[co_index+2] == 'C':  #  ...C(R)N(C
                            if s[co_index+3:co_index+5] != '=O':  #  ...C(R)N(C=O)
                                if s[co_index+3] == ')':  #  ...C(R)N(C)
                                    if s[co_index+4] == 'C':  #  ...C(R)N(C)C
                                        if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(C)C
                                            if s[co_index+5:co_index+7] != '=O':  #  ...C(R)N(C)C=O
                                                if len(s) > co_index+5 and s[co_index+5] == '(':  #  ...C(R)N(C)C(
                                                    if s[co_index+6:co_index+8] != '=O':  #  ...C(R)N(C)C(=O)
                                                        closing_parenthesis = find_closing_parenthesis(s, co_index+5)
                                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)N(C)C(...)=O
                                                            tertiary_amine_number = tertiary_amine_number + 1

                                                else:
                                                    tertiary_amine_number = tertiary_amine_number + 1


                                elif s[co_index+3] == '(':  #  ...C(R)N(C(
                                    if s[co_index+4:co_index+6] != '=O':  #  ...C(R)N(C(=O)
                                        closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)N(C(A)=O
                                            second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                            if s[second_closing_parenthesis+1] == 'C':  #  ...C(R)N(C(A)B)C
                                                if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(C(A)B)C
                                                    if s[second_closing_parenthesis+2:second_closing_parenthesis+4] != '=O':  #  ...C(R)N(C(A)B)C=O
                                                        if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+2] == '(':  #  ...CN(C(A)B)C(
                                                            if s[second_closing_parenthesis+3:second_closing_parenthesis+5] != '=O':  #  ...C(R)N(C(A)B)C(=O)
                                                                third_closing_parenthesis = find_closing_parenthesis(s, second_closing_parenthesis+2)
                                                                if s[third_closing_parenthesis+1:third_closing_parenthesis+3] != '=O':  #  ...C(R)N(C(A)B)C(I)=O
                                                                    tertiary_amine_number = tertiary_amine_number + 1

                                                        else:
                                                            tertiary_amine_number = tertiary_amine_number + 1


                                else:
                                    closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[closing_parenthesis+1] == 'C':  #  ...C(R)N(CA)C
                                        if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(R)N(CA)C
                                            if s[closing_parenthesis+2:closing_parenthesis+4] != '=O':  #  ...C(R)N(CA)C=O
                                                if len(s) > closing_parenthesis+2 and s[closing_parenthesis+2] == '(':  #  CN(CA)C(
                                                    if s[closing_parenthesis+3:closing_parenthesis+5] != '=O':  #  ...C(R)N(CA)C(=O)
                                                        second_closing_parenthesis = find_closing_parenthesis(s, closing_parenthesis+2)
                                                        if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O':  #  ...C(R)N(CA)C(...)=O
                                                            tertiary_amine_number = tertiary_amine_number + 1

                                                else:
                                                    tertiary_amine_number = tertiary_amine_number + 1


            elif s[first_opening_parenthesis-1] == ')':  #  ...)(R)N...
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  #  ...C(R)(R')N...

                    if len(s) > co_index+1 and s[co_index+1] == '(':  #  ...C(R)(R')N(
                        if s[co_index+2] == 'C':  #  ...C(R)(R')N(C
                            if s[co_index+3:co_index+5] != '=O':  #  ...C(R)(R')N(C=O)
                                if s[co_index+3] == ')':  #  ...C(R)(R')N(C)
                                    if s[co_index+4] == 'C':  #  ...C(R)(R')N(C)C
                                        if s[co_index+5:co_index+7] != '=O':  #  ...C(R)(R')N(C)C=O
                                            if len(s) > co_index+5 and s[co_index+5] == '(':  #  ...C(R)(R')N(C)C(
                                                if s[co_index+6:co_index+8] != '=O':  #  ...C(R)(R')N(C)C(=O)
                                                    closing_parenthesis = find_closing_parenthesis(s, co_index+5)
                                                    if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)(R')N(C)C(...)=O
                                                        tertiary_amine_number = tertiary_amine_number + 1

                                            else:
                                                tertiary_amine_number = tertiary_amine_number + 1


                                elif s[co_index+3] == '(':  #  ...C(R)(R')N(C(
                                    if s[co_index+4:co_index+6] != '=O':  #  ...C(R)(R')N(C(=O)
                                        closing_parenthesis = find_closing_parenthesis(s, co_index+3)
                                        if s[closing_parenthesis+1:closing_parenthesis+3] != '=O':  #  ...C(R)(R')N(C(A)=O
                                            second_closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                            if s[second_closing_parenthesis+1] == 'C':  #  ...C(R)(R')N(C(A)B)C
                                                if s[second_closing_parenthesis+2:second_closing_parenthesis+4] != '=O':  #  ...C(R)(R')N(C(A)B)C=O
                                                    if len(s) > second_closing_parenthesis+2 and s[second_closing_parenthesis+2] == '(':  #  ...CN(C(A)B)C(
                                                        if s[second_closing_parenthesis+3:second_closing_parenthesis+5] != '=O':  #  ...C(R)(R')N(C(A)B)C(=O)
                                                            third_closing_parenthesis = find_closing_parenthesis(s, second_closing_parenthesis+2)
                                                            if s[third_closing_parenthesis+1:third_closing_parenthesis+3] != '=O':  #  ...C(R)(R')N(C(A)B)C(I)=O
                                                                tertiary_amine_number = tertiary_amine_number + 1

                                                    else:
                                                        tertiary_amine_number = tertiary_amine_number + 1


                                else:
                                    closing_parenthesis = find_closing_parenthesis(s, co_index+1)
                                    if s[closing_parenthesis+1] == 'C':  #  ...C(R)(R')N(CA)C
                                        if s[closing_parenthesis+2:closing_parenthesis+4] != '=O':  #  ...C(R)(R')N(CA)C=O
                                            if len(s) > closing_parenthesis+2 and s[closing_parenthesis+2] == '(':  #  CN(CA)C(
                                                if s[closing_parenthesis+3:closing_parenthesis+5] != '=O':  #  ...C(R)(R')N(CA)C(=O)
                                                    second_closing_parenthesis = find_closing_parenthesis(s, closing_parenthesis+2)
                                                    if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O':  #  ...C(R)(R')N(CA)C(I)=O
                                                        tertiary_amine_number = tertiary_amine_number + 1

                                            else:
                                                tertiary_amine_number = tertiary_amine_number + 1


        co_index = s.find('N', co_index + 1)

    return tertiary_amine_number

###################   Aromatic_amine  #######################################################
def aromatic_amine_group(s):

    aromatic_amine_number = 0

    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_rings = []


    co_index = s.find('N', 0)
    while co_index != -1:
        if s[co_index-1] =='c' and co_index-1 != -1:  # ...cN...
            if s[co_index+1:co_index+4] !='(=O)':  # excluding nitro group...cN(=O)...
                aromatic_amine_number = aromatic_amine_number + 1

        elif s[co_index-2:co_index] in list_of_rings or s[co_index-3:co_index] in list_of_rings:  # ...c1N...   or ...c12N
            if s[co_index+1:co_index+5] !='(=O)':  # excluding nitro group...c1N(=O)...
                aromatic_amine_number = aromatic_amine_number + 1

        elif s[co_index-1] ==')':  #  )N
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'c':  # c(...)N
                if s[co_index+1:co_index+5] !='(=O)':  # excluding nitro group...c(...)N(=O)...
                    aromatic_amine_number = aromatic_amine_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings or s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  # c1(...)N
                if s[co_index+1:co_index+5] !='(=O)':  # excluding nitro group...c1(...)N(=O)...
                    aromatic_amine_number = aromatic_amine_number + 1


        elif s[co_index-1] =='(':  #  (N
            if s[co_index-2] == 'c':  #  c(N
                if s[co_index+1:co_index+5] !='(=O)':  # excluding nitro group...c(N(=O)...
                    aromatic_amine_number = aromatic_amine_number + 1

            elif s[co_index-3:co_index-1] in list_of_rings or s[co_index-4:co_index-1] in list_of_rings:  #  c1(N
                if s[co_index+1:co_index+5] !='(=O)':  # excluding nitro group...c1(N(=O)...
                    aromatic_amine_number = aromatic_amine_number + 1


        if len(s) >= co_index+2 and s[co_index+1] =='c':  # ...Nc...
            aromatic_amine_number = aromatic_amine_number + 1

        elif len(s) >= co_index+2 and s[co_index+1] =='(':  # ...N(...
            if s[co_index+2] =='c':  # ...N(c...
                aromatic_amine_number = aromatic_amine_number + 1

                index = s.find('(', co_index+1)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1] == 'c':  # ...N(c...)c
                    aromatic_amine_number = aromatic_amine_number + 1

            else:
                index = s.find('(', co_index+1)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1] == 'c':  # ...N(...)c
                    if s[co_index+2:co_index+5] !='=O)':  # excluding nitro group...N(=O)c...
                        aromatic_amine_number = aromatic_amine_number + 1


        co_index = s.find('N', co_index + 1)


    return aromatic_amine_number

###################   Primary_amide  #######################################################
def primary_amide_group(s):
    primary_amide_number = 0
    # primary amide as first characters of the SMILE:
    if s[0:6] == 'O=C(N)':
        primary_amide_number = primary_amide_number + 1

    if s[0:6] == 'NC(=O)':
        primary_amide_number = primary_amide_number + 1

    # primary amide as last characters of the SMILE:
    if s[-6:] == 'C(=O)N':
        primary_amide_number = primary_amide_number + 1

    elif s[-6:] == 'C(N)=O':
        primary_amide_number = primary_amide_number + 1

    if s[-1] == 'N':
        if s[-4:-1] == 'O=C':
            primary_amide_number = primary_amide_number + 1
        elif len(s) >= 2 and s[-2] == ')':
            index = s.rfind(')', 0, -1)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-3:opening_parenthesis] == 'O=C':
                primary_amide_number = primary_amide_number + 1

    if s[-2:] == '=O':
        if s[-4:-2] == 'NC' and len(s) == 4:  #  /NC=O/
            primary_amide_number = primary_amide_number + 1

        elif len(s) >= 3 and s[-3] == ')':
            index = s.rfind(')', 0, -2)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-2:opening_parenthesis] == 'NC'  and len(s) == index+3:  #  /NC(...)=O/:
                primary_amide_number = primary_amide_number + 1

    # primary amide as last characters of a branch in the SMILE:
    co_index = s.find('C(=O)N)', 0)
    while co_index != -1:
        primary_amide_number = primary_amide_number + 1
        co_index = s.find('C(=O)N)', co_index + 1)

    co_index = s.find('C(N)=O)', 0)
    while co_index != -1:
        primary_amide_number = primary_amide_number + 1
        co_index = s.find('C(N)=O)', co_index + 1)

    co_index = s.find('N)', 0)
    while co_index != -1:
        if s[co_index-3:co_index] == 'O=C':
            primary_amide_number = primary_amide_number + 1
        elif s[co_index-1] == ')':
            index = s.rfind(')', 0, co_index)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-3:opening_parenthesis] == 'O=C':
                primary_amide_number = primary_amide_number + 1

        co_index = s.find('N)', co_index + 1)

    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-2:co_index] == 'NC':
            primary_amide_number = primary_amide_number + 1
        elif s[co_index-1] == ')' and co_index-1 != -1:
            index = s.rfind(')', 0, co_index)
            opening_parenthesis = find_opening_parenthesis(s, index)
            if s[opening_parenthesis-2:opening_parenthesis] == 'NC':
                primary_amide_number = primary_amide_number + 1
        co_index = s.find('=O)', co_index + 1)

    return primary_amide_number

###################   Secondary_amide  #######################################################
def secondary_amide_group(s):

    secondary_amide_number = 0
    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings
    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # secondary amide as first characters of the SMILE:
    if s[0:3] == 'O=C':
        if len(s) > 3 and s[3] == 'N':  # O=CN
            if s[4:6] in ['C(','c(']:  # O=CNC(
                secondary_amide_number = secondary_amide_number + 1

        elif len(s) > 3 and s[3] == '1':  # O=C1
            if s[4:6] in ['NC','Nc']:  # O=C1NC
                secondary_amide_number = secondary_amide_number + 1

            elif s[-2:] == 'N1':  #  O=C1...N1
                if s[-3] == 'C':   # O=C1...CN1/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[-3] == ')':   #  O=C1...)N1/
                    index1 = s.rfind(')', 0, -2)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # O=C1...C(...)N1/
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # O=C1...)(...)N1/
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # O=C1...C(...)(...)N1/
                            secondary_amide_number = secondary_amide_number + 1


        elif len(s) > 3 and s[3] == '(':  # O=C(
            if s[4:6] in ['NC','Nc']:  # O=C(NC
                secondary_amide_number = secondary_amide_number + 1

            else:
                index = s.find('(', 3)
                closing_parenthesis = find_closing_parenthesis(s, index)
                if s[closing_parenthesis+1:closing_parenthesis+3] in ['NC','Nc']:  # O=C(...)NC
                    secondary_amide_number = secondary_amide_number + 1

    # secondary amide as last characters of the SMILE:
    if s[-2:] == '=O':
        if s[-4:-2] == 'NC':  #  NC=O/
            if len(s) >= 5 and s[-5] in ['C','c']:  #  CNC=O/
                secondary_amide_number = secondary_amide_number + 1

            elif s[-6:-4] in list_tot:  #  C1NC=O/
                secondary_amide_number = secondary_amide_number + 1

            elif s[-7:-4] in list_tot:  #  C12NC=O/
                secondary_amide_number = secondary_amide_number + 1

            elif len(s) >= 5 and s[-5] == ')':  #  )NC=O/
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)CNC=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)CNC=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)CNC=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-1] == ')':  #  )(...)CNC=O/
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)(...)CNC=O/
                        secondary_amide_number = secondary_amide_number + 1


        elif len(s) >= 3 and s[-3] == ')':  #  )=O/
            index = s.rfind(')', 0, -2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1:first_opening_parenthesis+3] in ['C(NC','C(Nc']:  #  C(NC...)=O/
                secondary_amide_number = secondary_amide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] == 'NC':  #  NC(...)=O/
                if s[first_opening_parenthesis-3] in ['C','c']:  #  CNC(...)=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1NC(...)=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12NC(...)=O/
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )NC(...)=O/
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)NC(...)=O/
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)NC(...)=O/
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)NC(...)=O/
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)NC(...)=O/
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] in ['C','c']:  #  C(...)(...)NC(...)=O/
                            secondary_amide_number = secondary_amide_number + 1


    # secondary amide as middle characters of the SMILE:
    co_index = s.find('NC(=O)', 0)
    while co_index != -1:
        if s[co_index-1] in ['C','c'] and co_index-1 != -1:  #  CNC(=O)
            secondary_amide_number = secondary_amide_number + 1

        elif s[co_index-2:co_index] in list_tot:  #  C1NC(=O)
            secondary_amide_number = secondary_amide_number + 1

        elif s[co_index-3:co_index] in list_tot:  #  C12NC(=O)
            secondary_amide_number = secondary_amide_number + 1

        elif s[co_index-1] == ')':  #  )NC(=O)
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)NC(=O)
                secondary_amide_number = secondary_amide_number + 1

            if s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)NC(=O)
                secondary_amide_number = secondary_amide_number + 1

            if s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)NC(=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  )(...)NC(=O)
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)(...)NC(=O)
                    secondary_amide_number = secondary_amide_number + 1

        co_index = s.find('NC(=O)', co_index + 1)


    co_index = s.find('C(=O)N', 0)
    while co_index != -1:
        if len(s) > co_index+6 and s[co_index+6] in ['C','c']:  #  C(=O)NC
            secondary_amide_number = secondary_amide_number + 1

        elif s[co_index+6:co_index+8] in list_tot:  #  C(=O)NC1
            secondary_amide_number = secondary_amide_number + 1

        co_index = s.find('C(=O)N', co_index + 1)


    co_index = s.find('(NC', 0)
    while co_index != -1:
        if s[co_index-1] in ['C','c']:  #  C(NC
            if s[co_index+3:co_index+6] == '=O)':  #  C(NC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3:co_index+7] == '(=O)':  #  C(NC(=O)...
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3] == '(':  #  C(NC(
                index = s.find('(', co_index+3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C(NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1


        elif s[co_index-2:co_index] in list_tot:  #  C1(NC
            if s[co_index+3:co_index+6] == '=O)':  #  C1(NC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3:co_index+7] == '(=O)':  #  C1(NC(=O)...
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3] == '(':  #  C1(NC(
                index = s.find('(', co_index+3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C1(NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1


        elif s[co_index-3:co_index] in list_tot:  #  C12(NC
            if s[co_index+3:co_index+6] == '=O)':  #  C12(NC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3:co_index+7] == '(=O)':  #  C12(NC(=O)...
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index+3] == '(':  #  C12(NC(
                index = s.find('(', co_index+3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C12(NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1


        elif s[co_index-1] == ')':  #  )(NC
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)(NC
                if s[co_index+3:co_index+6] == '=O)':  #  C(...)(NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3:co_index+7] == '(=O)':  #  C(...)(NC(=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3] == '(':  #  C(...)(NC(
                    index = s.find('(', co_index+3)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C(...)(NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)(NC
                if s[co_index+3:co_index+6] == '=O)':  #  C1(...)(NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3:co_index+7] == '(=O)':  #  C1(...)(NC(=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3] == '(':  #  C1(...)(NC(
                    index = s.find('(', co_index+3)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C1(...)(NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1


            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)(NC
                if s[co_index+3:co_index+6] == '=O)':  #  C12(...)(NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3:co_index+7] == '(=O)':  #  C12(...)(NC(=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[co_index+3] == '(':  #  C12(...)(NC(
                    index = s.find('(', co_index+3)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[first_closing_parenthesis+1:first_closing_parenthesis+4] == '=O)':  #  C12(...)(NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1

        co_index = s.find('(NC', co_index + 1)

    # secondary amide as last characters of a branch in the SMILE:
    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-2:co_index] == 'NC':  #  NC=O)
            if s[co_index-3] in ['C','c']:  #  CNC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index-4:co_index-2] in list_tot:  #  C1NC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index-5:co_index-2] in list_tot:  #  C12NC=O)
                secondary_amide_number = secondary_amide_number + 1

            elif s[co_index-3] == ')':  #  )NC=O)
                index = s.rfind(')', 0, co_index-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)NC=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-1] == ')':   #  )(...)NC=O)
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)(...)NC=O)
                        secondary_amide_number = secondary_amide_number + 1

        elif s[co_index-1] == ')':  #  )=O)
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis:first_opening_parenthesis+3] == '(NC':  #  (NC...)=O)
                if s[first_opening_parenthesis-1] in ['C','c']:  #  C(NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(NC...)=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(NC...)=O)
                    secondary_amide_number = secondary_amide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] == 'NC':  #  NC(...)=O)
                if s[first_opening_parenthesis-3] in ['C','c']:  #  CNC(...)=O
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12NC(...)=O)
                    secondary_amide_number = secondary_amide_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )NC(...)=O)
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)NC(...)=O)
                        secondary_amide_number = secondary_amide_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)NC(...)=O)
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] in ['C','c']:  #  C(...)(...)NC(...)=O)
                            secondary_amide_number = secondary_amide_number + 1

        co_index = s.find('=O)', co_index + 1)

    if cyc_number != None: # if there is any cycle
        for i in range(1, cyc_number+1):
            end_str = 'C(=O)N' + str(i)
            if s.endswith(end_str):
                index = s.find(str(i), 0)
                if s[index-1] == 'C':  # C1...C(=O)N1/
                    secondary_amide_number = secondary_amide_number + 1


        for i in range(1, cyc_number+1):
            end_str = 'C(N' + str(i) + ')=O'
            if s.endswith(end_str):
                index = s.find(str(i), 0)
                if s[index-1] == 'C':  # C1...C(N1)=O/
                    secondary_amide_number = secondary_amide_number + 1


    return secondary_amide_number

###################   Tertiary_amide  #######################################################
def tertiary_amide_group(s):
    tertiary_amide_number = 0

    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings
        list_N = ['N' + str(i) for i in range(1, cyc_number+1)]


    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []
        list_N = []

    # tertiary amide as first characters of the SMILE:
    if s[0:2] == 'O=':  #  O=
        if s[2:4] == 'C(':  #  O=C(
            if s[4:6] == 'N(':  #  O=C(N(
                if s[6] in ['C', 'c']:  #  O=C(N(C
                    index = s.find('(', 5)
                    closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[closing_parenthesis+1] in ['C', 'c']:  #  O=C(N(C...)C
                        tertiary_amide_number = tertiary_amide_number + 1

            elif s[4:6] in list_N:  #  O=C(N1
                if s[6] in ['C', 'c']:  #  O=C(N1C
                    N_index = s.find(s[5], 6)
                    if s[N_index-1] in ['C', 'c']:  #  O=C(N1C...C1
                        tertiary_amide_number = tertiary_amide_number + 1

            index = s.find('(', 3)
            closing_parenthesis = find_closing_parenthesis(s, index)
            if s[closing_parenthesis+1:closing_parenthesis+3] == 'N(':  #  O=C(...)N(
                if s[closing_parenthesis+3] in ['C', 'c']:  #  O=C(...)N(C
                    index2 = s.find('(', closing_parenthesis+2)
                    second_closing_parenthesis = find_closing_parenthesis(s, index2)
                    if s[second_closing_parenthesis+1] in ['C', 'c']:  #  O=C(...)N(C...)C
                        tertiary_amide_number = tertiary_amide_number + 1

            elif s[closing_parenthesis+1:closing_parenthesis+3] in list_N:  #  O=C(...)N1
                if len(s) > closing_parenthesis+3 and s[closing_parenthesis+3] in ['C', 'c']:  #  O=C(...)N1C
                    N_index = s.find(s[closing_parenthesis+2], closing_parenthesis+3)
                    if s[N_index-1] in ['C', 'c']:  #  O=C(...)N1C...C1
                        tertiary_amide_number = tertiary_amide_number + 1


        elif s[2:4] in list_of_nonarom_rings:  #  O=C1
            index = s.find(s[3], 4)
            if s[4:6] == 'N(':  #  O=C1N(
                if s[6] in ['C', 'c']:  #  O=C1N(C
                    index = s.find('(', 5)
                    closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[closing_parenthesis+1] in ['C', 'c']:  #  O=C1N(C...)C
                        tertiary_amide_number = tertiary_amide_number + 1

            elif s[4:6] in list_N:  #  O=C1N2
                if len(s) > 6 and s[6] in ['C', 'c']:  #  O=C1N2C
                    N_index = s.find(s[5], 6)
                    if s[N_index-1] in ['C', 'c']:  #  O=C1N2C...C2
                        tertiary_amide_number = tertiary_amide_number + 1

            if s[index-1] == 'N':  #  O=C1...N1
                if s[index-2] in ['C', 'c']:  #  O=C1...CN1
                    if len(s) > index+1 and s[index+1] in ['C', 'c']:  #  O=C1...CN1C
                        tertiary_amide_number = tertiary_amide_number + 1

                elif s[index-3:index-1] in list_tot:  #  O=C1...C2N1
                    if len(s) > index+1 and s[index+1] in ['C', 'c']:  #  O=C1...C2N1C
                        tertiary_amide_number = tertiary_amide_number + 1


        elif s[2:5] == 'CN(':  #  O=CN(
            if s[5] in ['C', 'c']:  #  O=CN(C
                index = s.find('(', 4)
                closing_parenthesis = find_closing_parenthesis(s, index)
                if s[closing_parenthesis+1] in ['C', 'c']:  #  O=CN(C...)C
                    tertiary_amide_number = tertiary_amide_number + 1

        elif s[2:5] == 'CN1':  #  O=CN1
            if s[5] in ['C', 'c']:  #  O=CN1C
                N_index = s.find(s[4], 5)
                if s[N_index-1] in ['C', 'c']:  #  O=CN1C...C1
                    tertiary_amide_number = tertiary_amide_number + 1


    # tertiary amide as last characters of the SMILE:
    if s[-2:] == '=O':  #  =O/
        if s[-3] == 'C':  #  C=O/
            if len(s) >= 4 and s[-4] == ')':  #  )C=O/
                index = s.rfind(')', 0, -3)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'N':  #  N(...)C=O/
                    if s[first_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C=O/
                        if s[first_opening_parenthesis-2] == 'C':  #  CN(C...)C=O/
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1N(C...)C=O/
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-2] == ')':  #  )N(C...)C=O/
                            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                            second_opening_parenthesis = find_opening_parenthesis(s, index1)
                            if s[second_opening_parenthesis-1] == 'C':  #  C(...)N(C...)C=O/
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)N(C...)C=O/
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-1] == ')':  #  (...)N(C...)C=O/
                                index2 = s.rfind(')', 0, second_opening_parenthesis)
                                third_opening_parenthesis = find_opening_parenthesis(s, index2)
                                if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C=O/
                                    tertiary_amide_number = tertiary_amide_number + 1


            elif s[-5:-3] in list_N:  #  N1C=O/
                if s[-6] in ['C', 'c']: #  CN1C=O/
                    N_index = s.find(s[-4], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C1...CN1C=O/
                        tertiary_amide_number = tertiary_amide_number + 1

                elif s[-7:-5] in list_tot: #  C2N1C=O/
                    N_index = s.find(s[-4], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C1...C2N1C=O/
                        tertiary_amide_number = tertiary_amide_number + 1

        elif s[-4:-2] in list_tot:  #  C1=O/
            if s[-5] == ')':  #  )C1=O/
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'N':  #  N(...)C1=O/
                    if s[first_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C1=O/
                        if s[first_opening_parenthesis-2] == 'C':  #  CN(C...)C1=O/
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C2N(C...)C1=O/
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-2] == ')':  #  )N(C...)C1=O/
                            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                            second_opening_parenthesis = find_opening_parenthesis(s, index1)
                            if s[second_opening_parenthesis-1] == 'C':  #  C(...)N(C...)C1=O/
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C2(...)N(C...)C1=O/
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-1] == ')':  #  (...)N(C...)C1=O/
                                index2 = s.rfind(')', 0, second_opening_parenthesis)
                                third_opening_parenthesis = find_opening_parenthesis(s, index2)
                                if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C1=O/
                                    tertiary_amide_number = tertiary_amide_number + 1


            elif s[-6:-4] in list_N:  #  N2C1=O/
                if s[-7] in ['C', 'c']: #  CN2C1=O/
                    N_index = s.find(s[-5], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C2...CN2C1=O/
                        tertiary_amide_number = tertiary_amide_number + 1

                elif s[-8:-6] in list_tot: #  C3N2C1=O/
                    N_index = s.find(s[-5], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C2...C3N2C1=O/
                        tertiary_amide_number = tertiary_amide_number + 1

            for i in range(1,cyc_number+1):
                index_N = s.find('N'+str(i), 0)
                if index_N != None:
                    if s[index_N+2] in ['C', 'c']:  #  N1C...C1=O
                        if s[index_N-1] in ['C', 'c']:  #  CN1C...C1=O
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[index_N-2:index_N] in list_tot:  #  C2N1C...C1=O
                            tertiary_amide_number = tertiary_amide_number + 1

    if s[-3:] == ')=O':  #  )=O
        index = s.rfind(')', 0, -2)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-1] == 'C':  #  C(...)=O
            if s[first_opening_parenthesis-2] == ')':  #  )C(...)=O
                index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                second_opening_parenthesis = find_opening_parenthesis(s, index1)
                if s[second_opening_parenthesis-1] == 'N':  #  N(...)C(...)=O
                    if s[second_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C(...)=O
                        if s[second_opening_parenthesis-2] == 'C':  #  CN(C...)C(...)=O
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[second_opening_parenthesis-3:second_opening_parenthesis-1] in list_tot:  #  C1N(C...)C(...)=O
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[second_opening_parenthesis-2] == ')':  #  )N(C...)C(...)=O
                            index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                            third_opening_parenthesis = find_opening_parenthesis(s, index2)
                            if s[third_opening_parenthesis-1] == 'C':  #  C(...)N(C...)C(...)=O
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[third_opening_parenthesis-2:third_opening_parenthesis] in list_tot:  #  C1(...)N(C...)C(...)=O
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[third_opening_parenthesis-1] == ')':  #  )(...)N(C...)C(...)=O
                                index3 = s.rfind(')', 0, third_opening_parenthesis)
                                forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                                if s[forth_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C(...)=O
                                    tertiary_amide_number = tertiary_amide_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_N:  #  N1C(...)=O/
                if s[first_opening_parenthesis-4] in ['C', 'c']: #  CN1C(...)=O/
                    N_index = s.find(s[first_opening_parenthesis-2], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C1...CN1C(...)=O/
                        tertiary_amide_number = tertiary_amide_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-3] in list_tot: #  C1N2C(...)=O/
                    N_index = s.find(s[first_opening_parenthesis-2], 0)
                    if s[N_index-1] in ['C', 'c']:  #  C1...C2N1C(...)=O/
                        tertiary_amide_number = tertiary_amide_number + 1


        elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] == 'C(N(':  #  C(N(...)...)=O
            if s[first_opening_parenthesis+3] in ['C', 'c']:  #  C(N(C...)...)=O
                index4 = s.find('(', first_opening_parenthesis+2)
                closing_parenthesis = find_closing_parenthesis(s, index4)
                if s[closing_parenthesis+1] in ['C', 'c']:  #  C(N(C...)...)=O
                    tertiary_amide_number = tertiary_amide_number + 1

    # tertiary amide as middle characters of the SMILE:
    co_index = s.find(')C(=O)', 0)
    while co_index != -1:
        index = s.rfind(')', 0, co_index+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-1] == 'N':  #  N(...)C(=O)
            if s[second_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C(=O)
                if s[second_opening_parenthesis-2] in ['C', 'c'] and second_opening_parenthesis-2 != -1: #  CN(C...)C(=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[second_opening_parenthesis-3:second_opening_parenthesis-1] in list_tot:  #  C1N(C...)C(=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[second_opening_parenthesis-2] == ')':  #  )N(C...)C(=O)
                    index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] in ['C', 'c']:  #  C(...)N(C...)C(=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[third_opening_parenthesis-2:third_opening_parenthesis] in list_tot:  #  C1(...)N(C...)C(=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[third_opening_parenthesis-1] == ')':  #  )(...)N(C...)C(=O)
                        index3 = s.rfind(')', 0, third_opening_parenthesis)
                        forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                        if s[forth_opening_parenthesis-1] == 'C':    #  C(...)(...)N(C...)C(=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                elif s[second_opening_parenthesis-2] == '(':  #  (N(C...)C(=O))
                    if s[second_opening_parenthesis-3] in ['C', 'c']:  #  C(N(C...)C(=O))
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-4:second_opening_parenthesis-2] in list_tot:  #  C1(N(C...)C(=O))
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-3] == ')':  #  )(N(C...)C(=O))
                        index2 = s.rfind(')', 0, second_opening_parenthesis-2)
                        third_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(N(C...)C(=O))
                            tertiary_amide_number = tertiary_amide_number + 1

        co_index = s.find(')C(=O)', co_index + 1)

    co_index1 = s.find('C(=O)N(', 0)
    while co_index1 != -1:
        if s[co_index1+7] in ['C', 'c']:  #  C(=O)N(C
            index = s.find('(', co_index1+6)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1] in ['C', 'c']:  #  C(=O)N(C...)C
                tertiary_amide_number = tertiary_amide_number + 1

        co_index1 = s.find('C(=O)N(', co_index1 + 1)

    co_index2 = s.find(')C=O)', 0)
    while co_index2 != -1:
        index = s.rfind(')', 0, co_index2+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-2:second_opening_parenthesis] == '(N':  #  (N(C...)C=O)
            if s[second_opening_parenthesis+1] in ['C', 'c']:  #  (N(C...)C=O)
                if s[second_opening_parenthesis-3] in ['C', 'c']:  #  C(N(C...)C=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[second_opening_parenthesis-4:second_opening_parenthesis-2] in list_tot:  #  C1(N(C...)C=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[second_opening_parenthesis-3] == ')':  #  )(N(C...)C=O)
                    index2 = s.rfind(')', 0, second_opening_parenthesis-2)
                    third_opening_parenthesis = find_opening_parenthesis(s, index2)
                    if s[third_opening_parenthesis-1] == 'C':  #  C(...)(N(C...)C=O)
                        tertiary_amide_number = tertiary_amide_number + 1

        co_index2 = s.find(')C=O)', co_index2 + 1)



    if cyc_number != None:
        for i in range(1,cyc_number+1):
            co_index2 = s.find(')C'+ str(i) +'=O)', 0)  # )C1=O)
            while co_index2 != -1:
                index = s.rfind(')', 0, co_index2+1)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-2:second_opening_parenthesis] == '(N':  #  (N(C...)C1=O)
                    if s[second_opening_parenthesis+1] in ['C', 'c']:  #  (N(C...)C1=O)
                        if s[second_opening_parenthesis-3] in ['C', 'c']:  #  C(N(C...)C1=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[second_opening_parenthesis-4:second_opening_parenthesis-2] in list_tot:  #  C2(N(C...)C1=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[second_opening_parenthesis-3] == ')':  #  )(N(C...)C1=O)
                            index2 = s.rfind(')', 0, second_opening_parenthesis-2)
                            third_opening_parenthesis = find_opening_parenthesis(s, index2)
                            if s[third_opening_parenthesis-1] == 'C':  #  C(...)(N(C...)C1=O)
                                tertiary_amide_number = tertiary_amide_number + 1


                co_index2 = s.find(')C'+ str(i) +'=O)', co_index2 + 1)


    co_index3 = s.find(')=O)', 0)
    while co_index3 != -1:
        index = s.rfind(')', 0, co_index3+1)
        second_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[second_opening_parenthesis-2:second_opening_parenthesis] == ')C':  #  )C(...)=O)
            index2 = s.rfind(')', 0, second_opening_parenthesis-1)
            third_opening_parenthesis = find_opening_parenthesis(s, index2)
            if s[third_opening_parenthesis-2:third_opening_parenthesis] == '(N':  #  (N(...)C(...)=O)
                if s[third_opening_parenthesis+1] in ['C', 'c']:  #  (N(C...)C(...)=O)
                    if s[third_opening_parenthesis-3] in ['C', 'c']:  #  C(N(C...)C(...)=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[third_opening_parenthesis-4:third_opening_parenthesis-2] in list_tot:  #  C1(N(C...)C(...)=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[third_opening_parenthesis-3] == ')':  #  )(N(C...)C(...)=O)
                        index3 = s.rfind(')', 0, third_opening_parenthesis-2)
                        forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                        if s[forth_opening_parenthesis-1] == 'C':  #  C(...)(N(C...)C(...)=O)
                            tertiary_amide_number = tertiary_amide_number + 1

        co_index3 = s.find(')=O)', co_index3 + 1)

    # tertiary amide as last characters of a branch in the SMILE:
    co_index4 = s.find(')C=O)', 0)
    while co_index4 != -1:
        index = s.rfind(')', 0, co_index4+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-1] == 'N':  #  N(...)C=O)
            if s[first_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C=O)
                if s[first_opening_parenthesis-2] in ['C', 'c']:  #  CN(C...)C=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1N(C...)C=O)
                    tertiary_amide_number = tertiary_amide_number + 1

                elif s[first_opening_parenthesis-2] == ')':  #  )N(C...)C=O)
                    index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                    second_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[second_opening_parenthesis-1] in ['C', 'c']:  #  C(...)N(C...)C=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)N(C...)C=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)N(C...)C=O)
                        index2 = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C=O)
                            tertiary_amide_number = tertiary_amide_number + 1

        co_index4 = s.find(')C=O)', co_index4 + 1)


    if cyc_number != None:
        for i in range(1,cyc_number+1):
            co_index4 = s.find(')C'+str(i)+'=O)', 0)  #  )C1=O)
            while co_index4 != -1:
                index = s.rfind(')', 0, co_index4+1)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'N':  #  N(...)C1=O)
                    if s[first_opening_parenthesis+1] in ['C', 'c']:  #  N(C...)C1=O)
                        if s[first_opening_parenthesis-2] in ['C', 'c']:  #  CN(C...)C1=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C2N(C...)C1=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[first_opening_parenthesis-2] == ')':  #  )N(C...)C1=O)
                            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
                            second_opening_parenthesis = find_opening_parenthesis(s, index1)
                            if s[second_opening_parenthesis-1] in ['C', 'c']:  #  C(...)N(C...)C1=O)
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C2(...)N(C...)C1=O)
                                tertiary_amide_number = tertiary_amide_number + 1

                            elif s[second_opening_parenthesis-1] == ')':  #  )(...)N(C...)C1=O)
                                index2 = s.rfind(')', 0, second_opening_parenthesis)
                                third_opening_parenthesis = find_opening_parenthesis(s, index2)
                                if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C1=O)
                                    tertiary_amide_number = tertiary_amide_number + 1

                co_index4 = s.find(')C'+str(i)+'=O)', co_index4 + 1)


    co_index5 = s.find(')=O)', 0)
    while co_index5 != -1:
        index = s.rfind(')', 0, co_index5+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-2:first_opening_parenthesis] == ')C':  #  )C(...)=O)
            index1 = s.rfind(')', 0, first_opening_parenthesis-1)
            second_opening_parenthesis = find_opening_parenthesis(s, index1)
            if s[second_opening_parenthesis-1] == 'N':   #  N(...)C(...)=O)
                if s[second_opening_parenthesis+1] in ['C', 'c']:   #  N(C...)C(...)=O)
                    if s[second_opening_parenthesis-2] in ['C', 'c']:  #  CN(C...)C(...)=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis-1] in list_tot:  #  C1N(C...)C(...)=O)
                        tertiary_amide_number = tertiary_amide_number + 1

                    elif s[second_opening_parenthesis-2] == ')':  #  )N(C...)C(...)=O)
                        index2 = s.rfind(')', 0, second_opening_parenthesis-1)
                        third_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[third_opening_parenthesis-1] in ['C', 'c']:  #  C(...)N(C...)C(...)=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[third_opening_parenthesis-2:third_opening_parenthesis] in list_tot:  #  C1(...)N(C...)C(...)=O)
                            tertiary_amide_number = tertiary_amide_number + 1

                        elif s[third_opening_parenthesis-1] == ')':  #  )(...)N(C...)C(...)=O)
                            index3 = s.rfind(')', 0, third_opening_parenthesis)
                            forth_opening_parenthesis = find_opening_parenthesis(s, index3)
                            if s[forth_opening_parenthesis-1] == 'C':  #  C(...)(...)N(C...)C(...)=O)
                                tertiary_amide_number = tertiary_amide_number + 1

        elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] == 'C(N(':  #  C(N(...)=O)
            if s[first_opening_parenthesis+3] in ['C', 'c']:  #  C(N(C...)=O)
                index4 = s.find('(', first_opening_parenthesis+3)
                closing_parenthesis = find_closing_parenthesis(s, index4)
                if s[closing_parenthesis+1] in ['C', 'c']:  #  C(N(C...)C...)=O)
                    tertiary_amide_number = tertiary_amide_number + 1

        co_index5 = s.find(')=O)', co_index5 + 1)


    return tertiary_amide_number

###################   Carbonylperoxynitrate  #######################################################
def carbonylperoxynitrate_group(s):
    carbonylperoxynitrate_number = 0
    # carbonylperoxynitrate as first characters of the SMILE
    if s[0:14] == '[O-][N+](=O)OO':  # /[O-][N+](=O)OO
        if s[14] == 'C':  # /[O-][N+](=O)OOC
            if s[15:17] == '=O':  # /[O-][N+](=O)OOC=O
                carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

            elif s[15] == '(':  # /[O-][N+](=O)OOC(
                if s[16:19] == '=O)':  # /[O-][N+](=O)OOC(=O)
                    carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

                else:
                    index = s.find('(', 15)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  # /[O-][N+](=O)OOC(...)=O/
                        carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1


    elif s[0:9] == 'O=N(=O)OO':  # /O=N(=O)OO
            if s[9] == 'C':  # /O=N(=O)OOC
                if s[10:12] == '=O':  # /O=N(=O)OOC=O
                    carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

                elif s[10] == '(':  # /O=N(=O)OOC(
                    if s[11:14] == '=O)':  # /O=N(=O)OOC(=O)
                        carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

                    else:
                        index = s.find('(', 10)
                        first_closing_parenthesis = find_closing_parenthesis(s, index)
                        if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O':  # /O=N(=O)OOC(...)=O/
                            carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1


    elif s[0:3] == 'O=C':  # O=C
        if s[3:17] == 'OO[N+](=O)[O-]':  #  /O=COO[N+](=O)[O-]
            carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

        elif s[3:12] == 'OON(=O)=O':  #  /O=COON(=O)=O/
            carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

        elif s[3] == '(':  # O=C(
            if s[4:18] == 'OO[N+](=O)[O-]':  # O=C(OO[N+](=O)[O-])...
                carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

            elif s[4:13] == 'OON(=O)=O':  # O=C(OON(=O)=O)...
                carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

            else:
                index = s.find('(', 3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+15] == 'OO[N+](=O)[O-]':  # O=C(...)OO[N+](=O)[O-]...
                    carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+10] == 'OON(=O)=O':  # O=C(...)OON(=O)=O...
                    carbonylperoxynitrate_number = carbonylperoxynitrate_number + 1


    # carbonylperoxynitrate as middle characters of the SMILE
    count_1 = s.count('C(=O)OON(=O)=O')
    count_2 = s.count('C(=O)OO[N+](=O)[O-]')
    count_3 = s.count('C(OON(=O)=O)=O')
    count_4 = s.count('C(OO[N+](=O)[O-])=O')

    carbonylperoxynitrate_number = carbonylperoxynitrate_number + count_1 + count_2 + count_3 + count_4


    return carbonylperoxynitrate_number
###################   Peroxide  #######################################################
def peroxide_group(s):

    peroxy_number = 0

    cyc_number = find_highest_digit(s)

    if cyc_number != None:
        list_of_rings = ['C' + str(i) for i in range(1, cyc_number+1)] + ['c' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_rings = []

    # Ether as middle characters of the SMILE
    co_index = s.find('OOC', 1)
    if co_index != 0:  # /OOC
        while co_index != -1:
            if s[co_index-1] in ['C','c']:  # COOC
                peroxy_number = peroxy_number + 1

            elif s[co_index-2:co_index] in list_of_rings:  # C1OOC
                peroxy_number = peroxy_number + 1

            elif s[co_index-3:co_index] in list_of_rings:  # C12OOC
                peroxy_number = peroxy_number + 1

            elif s[co_index-1] in ['(']:  # ...(OOC...
                if s[co_index-2] in ['C','c']:  # ...C(OOC...
                    peroxy_number = peroxy_number + 1

                elif s[co_index-3:co_index-1] in list_of_rings:  # ...C1(OOC...
                    peroxy_number = peroxy_number + 1

                elif s[co_index-4:co_index-1] in list_of_rings:  # ...C12(OOC...
                    peroxy_number = peroxy_number + 1

                elif s[co_index-2] == ')':  #  ...)(OOC...
                    index = s.rfind(')', 0, co_index-1)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)(OOC...
                        peroxy_number = peroxy_number + 1

                    elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  #  C1(...)(OOC...
                        peroxy_number = peroxy_number + 1

                    elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_of_rings:  #  C12(...)(OOC...
                        peroxy_number = peroxy_number + 1


            elif s[co_index-1] == ')':  #  ...)OOC...
                index = s.rfind(')', 0, co_index)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:  #  ...C(...)OOC...
                    peroxy_number = peroxy_number + 1

                elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_of_rings:  #  ...C1(...)OOC...
                    peroxy_number = peroxy_number + 1

                elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_of_rings:  #  ...C12(...)OOC...
                    peroxy_number = peroxy_number + 1

                elif s[second_opening_parenthesis-1] == ')':  #  ...)(...)OOC...
                    index = s.rfind(')', 0, second_opening_parenthesis)
                    third_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[third_opening_parenthesis-1] in ['C','c']:  #  ...C(...)(...)OOC...
                        peroxy_number = peroxy_number + 1

            co_index = s.find('OOC', co_index + 1)

    # APRAM SMILES   O(C...)OC...
    if s.startswith('O(C') or s.startswith('O(c'):
        first_closing_parenthesis = find_closing_parenthesis(s, 1)
        if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == 'OC' or s[first_closing_parenthesis+1:first_closing_parenthesis+3] == 'Oc':
            peroxy_number = peroxy_number + 1

    if cyc_number != None: # if there is any cycle
        # Ending with as a cycle
        for i in range(1, cyc_number+1):
            co_index = s.find('(OO' + str(i) + ')', 0)  #   ...(OO1)
            while co_index != -1:
                if s[co_index-1] == 'C':  # ...C(OO1)...
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...C(OO1)...
                        peroxy_number = peroxy_number + 1

                elif s[co_index-1] == ')':  # ...)(OO1)...
                    index1 = s.rfind(')', 0, co_index)
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)(OO1)...
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(OO1)...
                            peroxy_number = peroxy_number + 1

                co_index = s.find('(OO' + str(i) + ')', co_index + 1)

        for i in range(1, cyc_number+1):
            end_str = 'OO' + str(i)
            if s.endswith(end_str):
                if len(s) > len(end_str) and s[-len(end_str) - 1] == 'C':   # COO1/
                    index = s.find(str(i), 0)
                    if s[index-1] == 'C':  # C1...COO1/
                        peroxy_number = peroxy_number + 1

                elif len(s) > len(end_str) and s[-len(end_str) - 1] == ')':   # )OO1/
                    index1 = s.rfind(')', 0, -len(end_str))
                    first_opening_parenthesis = find_opening_parenthesis(s, index1)
                    if s[first_opening_parenthesis-1] == 'C':  # ...C(...)OO1/
                        index = s.find(str(i), 0)
                        if s[index-1] == 'C':  # C1...C(...)OO1/
                            peroxy_number = peroxy_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OO1/
                        index2 = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index2)
                        if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)OO1/
                            index = s.find(str(i), 0)
                            if s[index-1] == 'C':  # C1...C(...)(...)OO1/
                                peroxy_number = peroxy_number + 1


    return peroxy_number

###################   Hydroperoxide  #######################################################
def hydroperoxide_group(s):
    hydroperoxide_number = 0
    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings

    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # hydroperoxide as last characters of the SMILE:
    if s[-2:] == 'OO':  #  OO/
        if len(s) >= 3 and s[-3] == 'C':  #  COO/
            if s[-5:-3] != 'O=':  #  /O=COO/ !
                hydroperoxide_number = hydroperoxide_number + 1

        elif len(s) >= 3 and s[-3] == ')':  #  )OO/
            index = s.rfind(')', 0, -2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)OO/  or c(...c1...)OO/
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  #  C(=O)OO/ !
                    if s[first_opening_parenthesis-3:first_opening_parenthesis-1] != 'O=':  #  O=C(...)OO/ !
                        hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-1] == ')':  #  )(...)OO/
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OO/
                    hydroperoxide_number = hydroperoxide_number + 1

                elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)(...)OO/
                    hydroperoxide_number = hydroperoxide_number + 1

                elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)(...)OO/
                    hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)OO/  or c1(...)OO/
                hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)OO/  or c12(...)OO/
                hydroperoxide_number = hydroperoxide_number + 1

        elif s[-4:-2] in list_tot:  #  C1OO/  or c1OO/
            hydroperoxide_number = hydroperoxide_number + 1

        elif s[-5:-2] in list_tot:  #  C12OO/  or c12OO/
            hydroperoxide_number = hydroperoxide_number + 1


    #  hydroperoxide as first characters of the SMILE:
    #  excluding if the string starts with hydroperoxide i.e., /COO/ since they were included above
    if s.startswith('C('):  #  /C(
        if s[2:5] == 'OO)':  #  /C(OO)
            if s[5:7] == '=O':  #  /C(OO)=O!
                pass
            elif s[5] == '(':  #  /C(OO)(
                if s[6:8] == '=O':  #  /C(OO)(=O)!
                    pass
                elif s[6:9] == 'OO)':  #  /C(OO)(OO)
                    if s[9:11] != '=O':  #  /C(OO)(OO)=O !
                        hydroperoxide_number = hydroperoxide_number + 2

                else:
                    index = s.find('(', 5)
                    second_closing_parenthesis = find_closing_parenthesis(s, index)
                    if s[second_closing_parenthesis+1:second_closing_parenthesis+3] != '=O':  #  /C(OO)(...)=O!
                        hydroperoxide_number = hydroperoxide_number + 1

            else:
                hydroperoxide_number = hydroperoxide_number + 1

        elif s[2:4] == '=O':  #  /C(=O) !
            pass
        else:
            index = s.find('(', 1)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1:first_closing_parenthesis+5] == '(OO)':  #  /C(...)(OO)...
                if s[first_closing_parenthesis+5:first_closing_parenthesis+7] != '=O':  #  /C(...)(OO)=O/ !
                    hydroperoxide_number = hydroperoxide_number + 1

    if s.startswith('OOC')  or s.startswith('OOc'):  #  /OOC  or  #  /OOc1
        if s[3:5] == '=O':  #  /OOC=O!
            pass

        elif len(s) > 3 and s[3] == '(':  #  /OOC(
            if s[4:6] == '=O':  #  /OOC(=O)!
                pass

            else:
                index = s.find('(', 3)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+3] == '=O': #  /OOC(...)=O!
                    pass

                else:
                    hydroperoxide_number = hydroperoxide_number + 1  # /OOC(OO)... can have two hydroperoxide groups which is considered in next parts below.

        else:
            hydroperoxide_number = hydroperoxide_number + 1




           # nonaromatic ring-connected hydroperoxide

    for j in list_tot:
        if s.startswith(j + '('):  #  /C1(  or  /c1(
            if s[3:6] == 'OO)':  #  /C1(OO)
                if s[6] == '(':  #  /C1(OO)(
                    if s[7:10] == 'OO)':  #  /C1(OO)(OO)
                        hydroperoxide_number = hydroperoxide_number + 2

                    else:
                        hydroperoxide_number = hydroperoxide_number + 1

                else:
                    hydroperoxide_number = hydroperoxide_number + 1

            elif s[3:5] == '=O':  #  /C1(=O) !
                pass

            else:
                index = s.find('(', 2)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[first_closing_parenthesis+1:first_closing_parenthesis+5] == '(OO)':  #  /C1(...)(OO)...
                    hydroperoxide_number = hydroperoxide_number + 1

        if s.startswith('OO' + j):  #  /OOC1
                hydroperoxide_number = hydroperoxide_number  # This was coverd already above.


    #  hydroperoxide as middle characters of the SMILE:
    co_index = s.find('(OO)', 0)
    while co_index != -1:
        if s[co_index-1] in ['C','c']:  #  C(OO)
            if co_index-1 != 0:  #  /C(OO)!
                if s[co_index-3:co_index-1] == 'O=':  #  O=C(OO)!
                    pass

                elif s[co_index+4:co_index+6] == '=O':  #  ...C(OO)=O!
                    pass

                else:
                    hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index-2:co_index] in list_tot:  #  C1(OO)...
            if co_index-2 != 0:  #  /C1(OO)...!
                hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index-3:co_index] in list_tot:  #  C12(OO)...
            hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index-1] == ')':  #  )(OO)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] == 'C':  #  C(...)(OO)...
                if first_opening_parenthesis-1 != 0:  #  /C(...)(OO)...!
                    hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)(OO)...
                if first_opening_parenthesis-2 != 0:  #  /C1(...)(OO)...!
                    hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)(OO)...
                hydroperoxide_number = hydroperoxide_number + 1

        co_index = s.find('(OO)', co_index + 1)


    #  hydroperoxide as last characters of a branch in the SMILE:
    co_index2 = s.find('OO)', 0)
    while co_index2 != -1:
        if s[co_index2-1] == 'C':  # COO)
            hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index2-2:co_index2] in list_tot:  # C1OO)
            hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index2-3:co_index2] in list_tot:  # C12OO)
            hydroperoxide_number = hydroperoxide_number + 1

        elif s[co_index2-1] == ')':  # )OO)
            index = s.rfind(')', 0, co_index2)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-1] in ['C','c']:  # ...C(...)OO)  or ...c(...)OO)
                if s[first_opening_parenthesis+1:first_opening_parenthesis+3] != '=O':  # ...C(=O)OO)!
                    hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # ...C1(...)OO)
                hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # ...C12(...)OO)
                hydroperoxide_number = hydroperoxide_number + 1

            elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OO)
                index = s.rfind(')', 0, first_opening_parenthesis)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] == 'C':  # C(...)(...)OO)
                    hydroperoxide_number = hydroperoxide_number + 1

        co_index2 = s.find('OO)', co_index2 + 1)


    return hydroperoxide_number

###################   Peroxide acid  #######################################################
def peroxide_acid_group(s):
    peroxide_acid_number = 0
    # peroxide_acid as first characters of the SMILE:
    allowed_prefixes = ['O=C(OO)', 'OOC(=O)']
    for prefix in allowed_prefixes:
        if s.startswith(prefix):  #  /O=C(OO)... or /OOC(=O)...
            peroxide_acid_number = peroxide_acid_number + 1

    # peroxide_acid as last characters of the SMILE:
    if s[-7:] == 'C(=O)OO' or s[-7:] == 'C(OO)=O':  #  ...C(=O)OO/ or ...C(OO)=O/
        peroxide_acid_number = peroxide_acid_number + 1

    # peroxide_acid as middle characters of the SMILE:
    co_index = s.find('C(=O)OO)', 0)
    while co_index != -1:
        peroxide_acid_number = peroxide_acid_number + 1

        co_index = s.find('C(=O)OO)', co_index + 1)

    co_index = s.find('C(OO)=O)', 0)
    while co_index != -1:
        peroxide_acid_number = peroxide_acid_number + 1

        co_index = s.find('C(OO)=O)', co_index + 1)

    if s == 'O=COO' or s == 'OOC=O':  # O=COO
        peroxide_acid_number = peroxide_acid_number + 1

    return peroxide_acid_number

###################   Nitrophenol  #######################################################
def nitrophenol_group(s):

    nitrophenol_number = 0

    cyc_number = find_highest_digit(s)
    if cyc_number != None:
        list_of_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
    else:
        list_of_rings = []
    # hydroxyl as last characters of the SMILE
    if s[-1] == 'O':
        if s[-3:-1] in list_of_rings:  # ...c1O/
            # including nitro
            index1 = s.rfind(s[-3:-1], 0, -3)
            if 'c(N(=O)=O)' in s[index1:-3]:  # c1...c(N(=O)=O)...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif 'c([N+](=O)[O-])' in s[index1:-3]:  # c1...c(N(=O)=O)...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif s[index1-7:index1] == 'O=N(=O)':    # O=N(=O)c1...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif s[index1-10:index1] == 'O=[N+][O-]':    # O=[N+][O-]c1...c1O/
                nitrophenol_number = nitrophenol_number + 1

            elif s[index1-12:index1] == '[O-][N+](=O)':    # [O-][N+](=O)c1...c1O/
                nitrophenol_number = nitrophenol_number + 1


        elif s[-2] == ')':  # ...)O
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            found = False   # Flag to control breaking both loops
            if s[first_opening_parenthesis-1] == 'c':  # ...c(...)O
                for j in list_of_rings:
                    if found:  # If the condition was met, break the for loop
                        break
                    index1 = s.find(j, 0, first_opening_parenthesis-1)
                    index2 = s.find(j, first_opening_parenthesis)
                    if index2 == -1:  # there is another cycle: c1...c1...c()O
                        pass

                    else:
                        index3 = s.find('c(', index1, index2)
                        while index3 != -1:
                            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                            if s[index3+2:index3+10] == 'N(=O)=O)': # c1...c(..c(N(=O)=O)...c1...)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True  # Set flag to break both loops
                                break

                            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...c(..c([N+](=O)[O-])...c1...)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...c(...c(...c1...)N(=O)=O)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...c(...c(...c1...)[N+](=O)[O-])O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c(..c1...)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1([N+](=O)[O-])...c(..c1...)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index2+2:index2+9] == 'N(=O)=O':  # c1...c(..c1N(=O)=O)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1...c(..c1[N+](=O)[O-])O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            index3 = s.find('c(', index3+1, index2)


            elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_of_rings:  # c1(...)O
                index1 = s[first_opening_parenthesis-2:first_opening_parenthesis]
                index2 = s.find(s[first_opening_parenthesis-2:first_opening_parenthesis], first_opening_parenthesis)
                index3 = s.find('c(', index1, index2)
                while index3 != -1:
                    first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                    if s[index3+2:index3+10] == 'N(=O)=O)': # c1(..c(N(=O)=O)...c1...)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1(..c(N(=O)=O)...c1...)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1(...c(...c1...)N(=O)=O)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1(...c(...c1...)N(=O)=O)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    elif s[index2+2:index2+9] == 'N(=O)=O':  # c1(..c1N(=O)=O)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1(..c1N(=O)=O)O
                        nitrophenol_number = nitrophenol_number + 1
                        break

                    index3 = s.find('c(', index3+1, index2)


    #  hydroxyl as first characters of the SMILE
    if s.startswith('Oc1'):  # Oc1
        index2 = s.find('c1', 3)  # Oc1...c1
        index3 = s.find('c(', 2, index2)
        while index3 != -1:
            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
            if s[index3+2:index3+10] == 'N(=O)=O)': # Oc1...c(N(=O)=O)...c1...
                nitrophenol_number = nitrophenol_number + 1
                break

            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # Oc1...c([N+](=O)[O-])...c1...
                nitrophenol_number = nitrophenol_number + 1
                break

            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # Oc1...c(N(=O)=O...c1...)...
                nitrophenol_number = nitrophenol_number + 1
                break

            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # Oc1...c([N+](=O)[O-]...c1...)...
                nitrophenol_number = nitrophenol_number + 1
                break

            elif s[index2+2:index2+9] == 'N(=O)=O':  # Oc1...c1N(=O)=O
                nitrophenol_number = nitrophenol_number + 1
                break

            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # Oc1...c1[N+](=O)[O-]
                nitrophenol_number = nitrophenol_number + 1
                break

            index3 = s.find('c(', index3+1, index2)


    #  hydroxyl as middle characters of the SMILE
    co_index = s.find('c(O)', 0)
    while co_index != -1:  # ...c(O)...
        found = False
        for j in list_of_rings:
            if found:  # If the condition was met, break the for loop
                break
            index1 = s.find(j, 0, co_index)
            index2 = s.find(j, co_index)
            if index2 == -1:
                pass

            else:
                index3 = s.find('c(', index1, index2)
                while index3 != -1:
                    first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                    if s[index3+2:index3+10] == 'N(=O)=O)': # c1...c(N(=O)=O)...c(O)...c1...
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...c(N(=O)=O)...c(O)...c1...
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...c(...c(O)...c1...)N(=O)=O...
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...c(...c(O)...c1...)N(=O)=O...
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index2+2:index2+9] == 'N(=O)=O':  # c1...c(O)...c1N(=O)=O
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1...c(O)...c1N(=O)=O
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...c(O)...c1
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...c(O)...c1
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index1-7:index1] == 'O=N(=O)':  # O=N(=O)c1...c(O)...c1
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index1-10:index1] == 'O=[N+][O-]':  # O=N(=O)c1...c(O)...c1
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    elif s[index1-12:index1] == '[O-][N+](=O)':  # [O-][N+](=O)c1...c(O)...c1
                        nitrophenol_number = nitrophenol_number + 1
                        found = True
                        break

                    index3 = s.find('c(', index3+1, index2)

        co_index = s.find('c(O)', co_index + 1)

    for j in list_of_rings:
        co_index2 = s.find(j + '(O)', 0)
        while co_index2 != -1:  # c1(O)...
            index = s.find(j, co_index2+4)
            index3 = s.find('c(', co_index2, index)
            while index3 != -1:
                first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                if s[index3+2:index3+10] == 'N(=O)=O)': # c1(O)...c(N(=O)=O)...c1...
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1(O)...c([N+](=O)[O-])...c1...
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1(O)...c(...c1...)N(=O)=O...
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1(O)...c(...c1...)[N+](=O)[O-]...
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index2+2:index2+9] == 'N(=O)=O':  # c1(O)...c1N(=O)=O
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # c1(O)...c1N(=O)=O
                    nitrophenol_number = nitrophenol_number + 1
                    break

                index3 = s.find('c(', index3+1, index)


            co_index2 = s.find(j + '(O)', co_index2 + 1)

    #  hydroxyl as last characters of a branch in the SMILE
    co_index3 = s.find('O)', 0)
    while co_index3 != -1:
        if s[co_index3-2:co_index3] in list_of_rings:  # ...(...c1O)
            index1 = s.rfind(s[co_index3-2:co_index3], 0, co_index3-2)
            index3 = s.find('c(', index1, co_index3-2)
            while index3 != -1:
                first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                if s[index3+2:index3+10] == 'N(=O)=O)': # c1...(...c(N(=O)=O)...c1O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index3+2:index3+15] == '[N+](=O)[O-])': # c1...(...c(N(=O)=O)...c1O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # c1...(...c(...c1O...)N(=O)=O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # c1...(...c(...c1O...)N(=O)=O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index1+2:index1+11] == '(N(=O)=O)':  # c1(N(=O)=O)...(...c1O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                elif s[index1+2:index1+16] == '([N+](=O)[O-])':  # c1(N(=O)=O)...(...c1O)
                    nitrophenol_number = nitrophenol_number + 1
                    break

                index3 = s.find('c(', index3+1, co_index3-2)


        elif s[co_index3-1] == ')':  # ...)O)
            index = s.rfind(')', 0, co_index3)
            second_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[second_opening_parenthesis-1] == 'c':  # ...c(...)O)...
                found = False
                for j in list_of_rings:
                    if found:  # If the condition was met, break the for loop
                        break
                    index1 = s.find(j, 0, second_opening_parenthesis-1)
                    index2 = s.find(j, second_opening_parenthesis)
                    if index2 == -1:
                        pass

                    else:
                        index3 = s.find('c(', index1, index2)
                        while index3 != -1:
                            first_closing_parenthesis = find_closing_parenthesis(s, index3+1)
                            if s[index3+2:index3+10] == 'N(=O)=O)': # (...c1...c(N(=O)=O)...c(...c1...)O)...
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index3+2:index3+15] == '[N+](=O)[O-])': # (...c1...c([N+](=O)[O-])...c(...c1...)O)...
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+8] == 'N(=O)=O':  # (...c1...c(...c(...c1...)N(=O)=O)O)...
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[first_closing_parenthesis+1:first_closing_parenthesis+13] == '[N+](=O)[O-]':  # (...c1...c(...c(...c1...)N(=O)=O)O)...
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index2+2:index2+9] == 'N(=O)=O':  # (...c1...c(..c1N(=O)=O)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            elif s[index2+2:index2+14] == '[N+](=O)[O-]':  # (...c1...c(..c1N(=O)=O)O
                                nitrophenol_number = nitrophenol_number + 1
                                found = True
                                break

                            index3 = s.find('c(', index3+1, index2)


        co_index3 = s.find('O)', co_index3 + 1)


    return nitrophenol_number

###################   Nitroester  #######################################################
def nitroester_group(s):
    ester_number = 0
    nitroester_number = 0
    cyc_number = find_cycle_number(s)

    # Replace '[N+](=O)[O-]', 'O=[N+][O-]', and '[O-][N+](=O)' with 'N(=O)=O' and 'O=N(=O)' respectively:
    s = s.replace('[N+](=O)[O-]', 'N(=O)=O')
    s = s.replace('O=[N+][O-]', 'O=N(=O)')
    s = s.replace('[O-][N+](=O)', 'O=N(=O)')

    if cyc_number != None:
        list_of_nonarom_rings = ['C' + str(i) for i in range(1, cyc_number+1)]
        list_of_arom_rings = ['c' + str(i) for i in range(1, cyc_number+1)]
        list_tot = list_of_nonarom_rings + list_of_arom_rings
    else:
        list_of_nonarom_rings = []
        list_of_arom_rings = []
        list_tot = []

    # ester as first characters of the SMILE
    if s[0:5] == 'O=C(O':  # /O=C(O
        if s[5] in ['C','c']:  # /O=C(OC or # /O=C(Oc
            index = s.find('(', 3)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1] in ['C','c']:  # /O=C(OC...)C
                co_index = s.find('N(=O)=O', first_closing_parenthesis+1)
                if co_index != -1:  # /O=C(OC...)C...N(=O)=O...
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

    elif s[0:4] == 'O=C(':
        if s[4] in ['C','c']:  #  /O=C(C or /O=C(c
            index = s.find('(', 3)
            first_closing_parenthesis = find_closing_parenthesis(s, index)
            if s[first_closing_parenthesis+1:first_closing_parenthesis+3] in ['OC','Oc']:  #  /O=C(C...)OC
                co_index = s.find('N(=O)=O', 3, first_closing_parenthesis)
                if co_index != -1:  #  /O=C(C...N(=O)=O...)OC
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

    # Ester as last characters of the SMILE
    if s[-2:] == '=O':  #  )=O/
        if s[-3] == ')':
            index = s.rfind(')', 0, -1)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-2:first_opening_parenthesis+2] in ['OC(C', 'OC(c']:  #  OC(C...)=O/
                if s[first_opening_parenthesis-3] in ['C','c'] :  #  COC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  COC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1OC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  C1OC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12OC(C...)=O/
                    co_index = s.find('N(=O)=O', first_opening_parenthesis, index)  #  C12OC(C...N(=O)=O...)=O/
                    if co_index != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )OC(C...)=O
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(C...)=O
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)OC(C...)=O/
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C1(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)OC(C...)=O/
                        co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                        if co_index != -1:  #  C1(...)OC(C...N(=O)=O...)=O
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1


                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)OC(C...)=O/
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)OC(C...)=O/
                            co_index = s.find('N(=O)=O', first_opening_parenthesis, s.rfind(')', 0, -1))
                            if co_index != -1:
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

            elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] in ['C(OC','C(Oc']:  #  C(OC...)=O/
                if s[first_opening_parenthesis-2] in ['C','c']:  #  CC(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...CC(OC...)=O/ or O=N(=O)...CC(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1C(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...C1C(OC...)=O/ or O=N(=O)...C1C(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-1] in list_tot:  #  C12C(OC...)=O/
                    co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index != -1 or co_index1 != -1:  #  ...N(=O)=O...C12C(OC...)=O/ or O=N(=O)...C12C(OC...)=O/
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2] == ')':  #  )C(OC...)=O/
                    index = s.rfind(')', 0, first_opening_parenthesis-1)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C1(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C1(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)C(OC...)=O/
                        co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index != -1 or co_index1 != -1:  #  C12(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C12(...)C(OC...)=O/
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)C(OC...)=O/
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)C(OC...)=O/
                            co_index = s.find('N(=O)=O', 0, first_opening_parenthesis)
                            co_index1 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                            if co_index != -1 or co_index1 != -1:  #  C(...)(...N(=O)=O...)C(OC...)=O/ or O=N(=O)...C(...)(...)C(OC...)=O/
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1


        elif s[-4:-2] == 'OC':  # ...OC=O
            if len(s) >= 5 and s[-5] == 'C':  # ...COC=O
                ester_number = ester_number + 1

            elif s[-6:-4] in list_tot:  # ...C1OC=O  or  # ...c1OC=O
                ester_number = ester_number + 1

            elif s[-7:-4] in list_tot:  # ...C12OC=O  or  # ...c12OC=O
                ester_number = ester_number + 1

            elif len(s) >= 5 and s[-5] == ')':  # ...)OC=O
                index = s.rfind(')', 0, -4)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  # ...C(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # ...C1(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # ...C12(...)OC=O
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OC=O
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)OC=O
                        ester_number = ester_number + 1


        elif s[-4:-2] in list_tot:  # ...C1=O
            if s[-5] == 'O': # ...OC1=O
                if s[-6] == 'C':  # ...COC1=O
                    ring_index = s.find(s[-3], 0, -5)
                    if s[ring_index-1] in ['C','c']:  #  C1...COC1=O/
                        ester_number = ester_number + 1

                elif s[-6] == ')':  # ...)OC1=O
                    index = s.rfind(')', 0, -5)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC1=O/
                        ring_index = s.find(s[-3], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)OC1=O/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC1=O/
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC1=O/
                            ring_index = s.find(s[-3], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)OC1=O/
                                ester_number = ester_number + 1



        elif s[-5:-2] in list_tot:  # ...C12=O
            if s[-6] == 'O': # ...OC12=O
                if s[-7] == 'C':  # ...COC21=O
                    ring_index = s.find(s[-4:-2], 0, -5)
                    if s[ring_index-1] in ['C','c']:  #  C12...COC12=O/
                        ester_number = ester_number + 1

                elif s[-7] == ')':  # ...)OC12=O
                    index = s.rfind(')', 0, -5)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC12=O/
                        ring_index = s.find(s[-4:-2], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)OC1=O/
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC12=O/
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC12=O/
                            ring_index = s.find(s[-4:-2], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)OC12=O/
                                ester_number = ester_number + 1



    # Ester as last characters of a channel in the SMILE
    co_index = s.find('=O)', 0)
    while co_index != -1:
        if s[co_index-1] == ')':  # ...)=O)...
            index = s.rfind(')', 0, co_index)
            first_opening_parenthesis = find_opening_parenthesis(s, index)
            if s[first_opening_parenthesis-2:first_opening_parenthesis+2] in ['OC(C', 'OC(c']:  #  OC(C...)=O)
                if s[first_opening_parenthesis-3] in ['C','c'] :  #  COC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  COC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  #  C1OC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  C1OC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  #  C12OC(C...)=O)
                    co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                    if co_index1 != -1:  #  C12OC(C...N(=O)=O...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3] == ')':  #  )OC(C...)=O)
                    index = s.rfind(')', 0, first_opening_parenthesis-2)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C1(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)OC(C...)=O)
                        co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                        if co_index1 != -1:  #  C1(...)OC(C...N(=O)=O...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)OC(C...)=O)
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)OC(C...)=O)
                            co_index1 = s.find('N(=O)=O', first_opening_parenthesis, index)
                            if co_index1 != -1:  #  C(...)(...)OC(C...N(=O)=O...)=O)
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

            elif s[first_opening_parenthesis-1:first_opening_parenthesis+3] in ['C(OC','C(Oc']:  #  C(OC...)=O)
                if s[first_opening_parenthesis-2] in ['C','c']:  #  CC(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...CC(OC...)=O) or ...N(=O)=O...CC(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis-1] in list_tot:  #  C1C(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C1C(OC...)=O) or ...N(=O)=O...C1C(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-4:first_opening_parenthesis-1] in list_tot:  #  C12C(OC...)=O)
                    co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C12C(OC...)=O) or ...N(=O)=O...C12C(OC...)=O)
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2] == ')':  #  )C(OC...)=O)
                    index = s.rfind(')', 0, first_opening_parenthesis-1)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C(...)C(OC...)=O) or ...N(=O)=O...C(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-2:second_opening_parenthesis] in list_tot:  #  C1(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C1(...)C(OC...)=O) or ...N(=O)=O...C1(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1


                    elif s[second_opening_parenthesis-3:second_opening_parenthesis] in list_tot:  #  C12(...)C(OC...)=O)
                        co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                        co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                        if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C12(...)C(OC...)=O) or ...N(=O)=O...C1(...)C(OC...)=O)
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

                    elif s[second_opening_parenthesis-1] == ')':  #  )(...)C(OC...)=O)
                        index = s.rfind(')', 0, second_opening_parenthesis)
                        third_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[third_opening_parenthesis-1] == 'C':  #  C(...)(...)C(OC...)=O)
                            co_index1 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                            co_index2 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                            if co_index1 != -1 or co_index2 != -1:  #  O=N(=O)...C(...)(...)C(OC...)=O) or ...N(=O)=O...C(...)(...)C(OC...)=O)
                                nitroester_number = nitroester_number + 1
                            else:
                                ester_number = ester_number + 1

        elif s[co_index-2:co_index] == 'OC':  # ...OC=O)
            if s[co_index-3] == 'C':  # ...COC=O
                ester_number = ester_number + 1

            elif s[co_index-4:co_index-2] in list_tot:  # ...C1OC=O)  or  # ...c1OC=O)
                ester_number = ester_number + 1

            elif s[co_index-5:co_index-2] in list_tot:  # ...C12OC=O)  or  # ...c12OC=O)
                ester_number = ester_number + 1

            elif s[co_index-3] == ')':  # ...)OC=O)
                index = s.rfind(')', 0, co_index-2)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] == 'C':  # ...C(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # ...C1(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # ...C12(...)OC=O)
                    ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] == ')':  # ...)(...)OC=O)
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] == 'C':  # ...C(...)(...)OC=O)
                        ester_number = ester_number + 1


        elif s[co_index-2:co_index] in list_tot:  # ...C1=O)
            if s[co_index-3] == 'O': # ...OC1=O)
                if s[co_index-4] == 'C':  # ...COC1=O)
                    ring_index = s.find(s[co_index-1], 0, co_index-4)
                    if s[ring_index-1] in ['C','c']:  #  C1...COC1=O)
                        ester_number = ester_number + 1

                elif s[co_index-4] == ')':  # ...)OC1=O)
                    index = s.rfind(')', 0, co_index-3)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC1=O)
                        ring_index = s.find(s[co_index-1], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C1...C(...)OC1=O)
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC1=O)
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC1=O)
                            ring_index = s.find(s[co_index-1], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C1...C(...)(...)OC1=O)
                                ester_number = ester_number + 1


        elif s[co_index-3:co_index] in list_tot:  # ...C12=O)
            if s[co_index-4] == 'O': # ...OC12=O)
                if s[co_index-5] == 'C':  # ...COC12=O)
                    ring_index = s.find(s[co_index-2:co_index], 0, co_index-4)
                    if s[ring_index-1] in ['C','c']:  #  C12...COC12=O)
                        ester_number = ester_number + 1

                elif s[co_index-5] == ')':  # ...)OC12=O)
                    index = s.rfind(')', 0, co_index-3)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)OC12=O)
                        ring_index = s.find(s[co_index-2:co_index], 0, first_opening_parenthesis)
                        if s[ring_index-1] in ['C','c']:  #  C12...C(...)OC12=O)
                            ester_number = ester_number + 1

                    elif s[first_opening_parenthesis-1] == ')':  #  )(...)OC12=O)
                        index = s.rfind(')', 0, first_opening_parenthesis)
                        second_opening_parenthesis = find_opening_parenthesis(s, index)
                        if s[second_opening_parenthesis-1] == 'C':  #  C(...)(...)OC12=O)
                            ring_index = s.find(s[co_index-2:co_index], 0, second_opening_parenthesis)
                            if s[ring_index-1] in ['C','c']:  #  C12...C(...)(...)OC12=O)
                                ester_number = ester_number + 1



        co_index = s.find('=O)', co_index + 1)

    # Ester as middle characters in the SMILE:
    co_index1 = s.find('C(=O)O', 0)
    while co_index1 != -1:
        if len(s) > co_index1+6 and s[co_index1+6] in ['C','c']:  # C(=O)OC
            if s[co_index1-1] in ['C','c']:  # CC(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-2:co_index1] in list_tot:  # C1C(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-3:co_index1] in list_tot:  # C12C(=O)OC
                co_index2 = s.find('N(=O)=O', 0, co_index1)
                co_index3 = s.find('O=N(=O)', 0, co_index1)
                if co_index2 != -1 or co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index1-1] == ')':  # )C(=O)OC
                index = s.rfind(')', 0, co_index1)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] in ['C', 'c']:  # C(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  # C1(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  # C12(...)C(=O)OC
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    if co_index2 != -1 or co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] in [')']:  # )(...)C(=O)OC
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C','c']: # C(...)(...)C(=O)OC
                        co_index2 = s.find('N(=O)=O', 0, co_index1)
                        co_index3 = s.find('O=N(=O)', 0, co_index1)
                        if co_index2 != -1 or co_index3 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

            elif s[co_index1-1] == '(':  # (C(=O)OC
                if s[co_index1-2] in ['C','c']:  # C(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-3:co_index1-1] in list_tot:  # C1(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-4:co_index1-1] in list_tot:  # C12(C(=O)OC
                    index = s.find('(', co_index1-1)
                    first_closing_parenthesis = find_closing_parenthesis(s, index)
                    co_index2 = s.find('N(=O)=O', 0, co_index1)
                    co_index3 = s.find('O=N(=O)', 0, co_index1)
                    co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                    if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index1-2] == ')':  # )(C(=O)OC
                    index = s.rfind(')', 0, co_index1-1)
                    third_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[third_opening_parenthesis-1] == 'C':  # C(...)(C(=O)OC
                        index2 = s.find('(', co_index1-1)
                        first_closing_parenthesis = find_closing_parenthesis(s, index2)
                        co_index2 = s.find('N(=O)=O', 0, co_index1)
                        co_index3 = s.find('O=N(=O)', 0, co_index1)
                        co_index4 = s.find('N(=O)=O', first_closing_parenthesis)
                        if co_index2 != -1 or co_index3 != -1 or co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

        co_index1 = s.find('C(=O)O', co_index1 + 1)

    co_index2 = s.find(')=O)', 0)
    while co_index2 != -1:
        index = s.rfind(')', 0, co_index2+1)
        first_opening_parenthesis = find_opening_parenthesis(s, index)
        if s[first_opening_parenthesis-2:first_opening_parenthesis+3] in ['(C(OC','(C(Oc']:  # (C(OC...)=O)
            if s[first_opening_parenthesis-3] in ['C','c']:  # C(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-4:first_opening_parenthesis-2] in list_tot:  # C1(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-5:first_opening_parenthesis-2] in list_tot:  # C12(C(OC...)=O)
                co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                co_index5 = s.find('N(=O)=O', co_index2)
                if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-3] == ')': # )(C(OC...)=O)
                index = s.rfind(')', 0, first_opening_parenthesis-2)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:   # C(...)(C(OC...)=O)
                    co_index3 = s.find('N(=O)=O', 0, first_opening_parenthesis)
                    co_index4 = s.find('O=N(=O)', 0, first_opening_parenthesis)
                    co_index5 = s.find('N(=O)=O', co_index2)
                    if co_index3 != -1 or co_index4 != -1 or co_index5 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

        elif s[first_opening_parenthesis-3:first_opening_parenthesis+2] in ['(OC(C','(OC(c']:
            if s[first_opening_parenthesis-4] in ['C','c']:  #  C(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-5:first_opening_parenthesis-3] in list_tot:  #  C1(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-6:first_opening_parenthesis-3] in list_tot:  #  C12(OC(C...)=O)
                co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                if co_index3 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[first_opening_parenthesis-4] == ')':  #  )(OC(C...)=O)
                index = s.rfind(')', 0, first_opening_parenthesis-3)
                second_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[second_opening_parenthesis-1] in ['C','c']:  #  C(...)(OC(C...)=O)
                    co_index3 = s.find('N(=O)=O', first_opening_parenthesis, co_index2)
                    if co_index3 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

        co_index2 = s.find(')=O)', co_index2 + 1)

    co_index3 = s.find('OC(=O)', 0)
    while co_index3 != -1:
        if s[co_index3+6] in ['C','c']:  #  OC(=O)C...
            if s[co_index3-1] in ['C','c']:  #  COC(=O)C...
                co_index4 = s.find('N(=O)=O', co_index3)
                if co_index4 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index3-2:co_index3] in list_tot:  #  C1OC(=O)C...
                co_index4 = s.find('N(=O)=O', co_index3)
                if co_index4 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index3-3:co_index3] in list_tot:  #  C12OC(=O)C...
                co_index4 = s.find('N(=O)=O', co_index3)
                if co_index4 != -1:
                    nitroester_number = nitroester_number + 1
                else:
                    ester_number = ester_number + 1

            elif s[co_index3-1] == ')':  #  )OC(=O)C...
                index = s.rfind(')', 0, co_index3)
                first_opening_parenthesis = find_opening_parenthesis(s, index)
                if s[first_opening_parenthesis-1] in ['C','c']:  #  C(...)OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-2:first_opening_parenthesis] in list_tot:  #  C1(...)OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-3:first_opening_parenthesis] in list_tot:  #  C12(...)OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[first_opening_parenthesis-1] in [')']:  #  )(...)OC(=O)C...
                    index = s.rfind(')', 0, first_opening_parenthesis)
                    second_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[second_opening_parenthesis-1] in ['C']:  #  C(...)(...)OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

            elif s[co_index3-1] == '(':  #  (OC(=O)C...
                index = s.find('(', co_index3-1)
                first_closing_parenthesis = find_closing_parenthesis(s, index)
                if s[co_index3-2] in ['C','c']:  #  C(OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-3:co_index3-1] in list_tot:  #  C1(OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-4:co_index3-1] in list_tot:  #  C12(OC(=O)C...
                    co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                    if co_index4 != -1:
                        nitroester_number = nitroester_number + 1
                    else:
                        ester_number = ester_number + 1

                elif s[co_index3-2] == ')':  #  )(OC(=O)C...
                    index = s.rfind(')', 0, co_index3-1)
                    first_opening_parenthesis = find_opening_parenthesis(s, index)
                    if s[first_opening_parenthesis-1] == 'C':  #  C(...)(OC(=O)C...
                        co_index4 = s.find('N(=O)=O', co_index3, first_closing_parenthesis)
                        if co_index4 != -1:
                            nitroester_number = nitroester_number + 1
                        else:
                            ester_number = ester_number + 1

        co_index3 = s.find('OC(=O)', co_index3 + 1)

    return nitroester_number

##########################################################################
##########################################################################
##########################################################################

def simpol_saturation_pressure(filename_in, filename_out,Antoine=True,T=298.15):
    for _ in simpol_saturation_pressure_gen(filename_in, filename_out,Antoine=Antoine,T=T,yieldABC=False):
        pass

def simpol_saturation_pressure_gen(filename_in, filename_out,Antoine=True,T=298.15,yieldABC=False):

    # Define B matrix directly in the code
    B = [
        [-4.26938E+02, 2.89223E-01, 4.42057E-03, 2.92846E-01],
        [-4.11248E+02, 8.96919E-01, -2.48607E-03, 1.40312E-01],
        [-1.46442E+02, 1.54528E+00, 1.71021E-03, -2.78291E-01],
        [3.50262E+01, -9.20839E-01, 2.24399E-03, -9.36300E-02],
        [-8.72770E+01, 1.78059E+00, -3.07187E-03, -1.04341E-01],
        [5.73335E+00, 1.69764E-02, -6.28957E-04, 7.55434E-03],
        [-2.61268E+02, -7.63282E-01, -1.68213E-03, 2.89038E-01],
        [-7.25373E+02, 8.26326E-01, 2.50957E-03, -2.32304E-01],
        [-7.29501E+02, 9.86017E-01, -2.92664E-03, 1.78077E-01],
        [-1.37456E+01, 5.23486E-01, 5.50298E-04, -2.76950E-01],
        [-7.98796E+02, -1.09436E+00, 5.24132E-03, -2.28040E-01],
        [-3.93345E+02, -9.51778E-01, -2.19071E-03, 3.05843E-01],
        [-1.44334E+02, -1.85617E+00, -2.37491E-05, 2.88290E-01],
        [4.05265E+01, -2.43780E+00, 3.60133E-03, 9.86422E-02],
        [-7.07406E+01, -1.06674E+00, 3.73104E-03, -1.44003E-01],
        [-7.83648E+02, -1.03439E+00, -1.07148E-03, 3.15535E-01],
        [-5.63872E+02, -7.18416E-01, 2.63016E-03, -4.99470E-02],
        [-4.53961E+02, -3.26105E-01, -1.39780E-04, -3.93916E-02],
        [3.71375E+01, -2.66753E+00, 1.01483E-03, 2.14233E-01],
        [-5.03710E+02, 1.04092E+00, -4.12746E-03, 1.82790E-01],
        [-3.59763E+01, -4.08458E-01, 1.67264E-03, -9.98919E-02],
        [-6.09432E+02, 1.50436E+00, -9.09024E-04, -1.35495E-01],
        [-1.02367E+02, -7.16253E-01, -2.90670E-04, -5.88556E-01],
        [-1.93802E+03, 6.48262E-01, 1.73245E-03, 3.47940E-02],
        [-5.26919E+00, 3.06435E-01, 3.25397E-03, -6.81506E-01],
        [-2.84042E+02, -6.25424E-01, -8.22474E-04, -8.80549E-02],
        [1.50093E+02, 2.39875E-02, -3.37969E-03, 1.52789E-02],
        [-2.03387E+01, -5.48718E+00, 8.39075E-03, 1.07884E-01],
        [-8.38064E+02, -1.09600E+00, -4.24385E-04, 2.81812E-01],
        [-5.27934E+01, -4.63689E-01, -5.11647E-03, 3.84965E-01],
        [-1.61520E+03, 9.01669E-01, 1.44536E-03, 2.66889E-01]
    ]

    file_out = open(filename_out, 'w')
    # Headline of output file:
    output_line = f"P_sat (Pa) at {T}K, Delta_H (kJ/mol) at {T}K, SMILES, carbon_number, ASA_carbon_number, aromatic_ring, non_aromatic_ring," \
                        f"double_bound_nonaromatic_ring, nonaromatic_CCCO, hydroxyl, aldehyde, ketone, carboxylic_acid," \
                        f"ester, ether, alicyclic_ether, aromatic_ether, nitrate, nitro, aromatic_hydroxyl, primary_amine," \
                        f"secondary_amine, tertiary_amine, aromatic_amine, primary_amide, secondary_amide, tertiary_amide," \
                        f"carbonylperoxynitrate, peroxide, hydroperoxide, peroxide_acid, nitrophenol, nitroester," \
                        f"log P_sat(T) [atm], Antoine_Eq for logP_sat(T) [atm], Enthalpy of vaporization (kJ/mol)\n"

    file_out.write(output_line)

    inorg_list = ['[C-]#[O+]', '[HH]', 'OO', '[N+](=O)(O)[O-]', 'O[O]', '[N+](=O)([O-])OO', 'N(=O)O',
              'OS(=O)[O]', '[N+](=O)([O-])O[N+](=O)[O-]', '[N]=O', 'N(=O)[O]', '[N+](=O)([O-])[O]',
              '[O]', '[O]', '[O-][O+]=O', '[OH]', 'O=S=O', 'O=S(=O)=O', 'nan']

    s_number = 1
    if isinstance(filename_in, str):
        file_in = open(filename_in, 'r')
        lines = file_in.readlines()
    else:
        lines = list(np.ravel([filename_in]))
    for line in lines:
        line = line.strip()  # Remove any leading or trailing whitespace
        line_intact = line
        line = line.replace('/', '').replace('\\', '').replace('[C]', 'C')
        line = line.strip()  # Remove any leading or trailing whitespace
        print('line = ', s_number)
        print('SMILES: ', line)
        column2  = carbon_number(line)
        column7  = ASA_carbon_number(line)
        column8  = aromatic_ring(line)
        column9  = non_aromatic_ring(line)
        column10 = double_bound_nonaromatic_carbons(line)
        column11 = nonaromatic_CCCO(line)
        column12 = hydroxyl_group(line)
        column13 = aldehyde_group(line)
        column14 = ketone_group(line)
        column15 = carboxylic_acid_group(line)
        column16 = ester_group(line)
        column17 = ether_group(line)
        column18 = alicyclic_ether(line)
        column19 = aromatic_ether_group(line)
        column20 = nitrate_number(line)
        column21 = nitro_group(line)
        column22 = aromatic_hydroxyl_group(line)
        column23 = primary_amine_group(line)
        column24 = secondary_amine_group(line)
        column25 = tertiary_amine_group(line)
        column26 = aromatic_amine_group(line)
        column27 = primary_amide_group(line)
        column28 = secondary_amide_group(line)
        column29 = tertiary_amide_group(line)
        column30 = carbonylperoxynitrate_group(line)
        column31 = peroxide_group(line)
        column32 = hydroperoxide_group(line)
        column33 = peroxide_acid_group(line)
        column34 = nitrophenol_group(line)
        column35 = nitroester_group(line)

        v = [1,
        carbon_number(line),
        ASA_carbon_number(line),
        aromatic_ring(line),
        non_aromatic_ring(line),
        double_bound_nonaromatic_carbons(line),
        nonaromatic_CCCO(line),
        hydroxyl_group(line),
        aldehyde_group(line),
        ketone_group(line),
        carboxylic_acid_group(line),
        ester_group(line),
        ether_group(line),
        alicyclic_ether(line),
        aromatic_ether_group(line),
        nitrate_number(line),
        nitro_group(line),
        aromatic_hydroxyl_group(line),
        primary_amine_group(line),
        secondary_amine_group(line),
        tertiary_amine_group(line),
        aromatic_amine_group(line),
        primary_amide_group(line),
        secondary_amide_group(line),
        tertiary_amide_group(line),
        carbonylperoxynitrate_group(line),
        peroxide_group(line),
        hydroperoxide_group(line),
        peroxide_acid_group(line),
        nitrophenol_group(line),
        nitroester_group(line)
        ]

        fun_num = sum(v[3:])

        b = [row[0]/T + row[1] + row[2]*T + row[3]*np.log(T) for row in B]

        h = [-2.303 * 8.314462618 * 0.001 * (row[0] - row[2] * (T**2) - row[3] * T) for row in B]

        result = []
        result_prime = []

        enthalpy = []
        enthalpy_prime = []

        for i in range(len(v)):
            result.append(v[i] * b[i])
            result_prime.append([v[i] * val for val in B[i]])
            enthalpy.append(v[i] * h[i])
            enthalpy_prime.append([v[i] * (-2.303 * 8.314462618) * val for val in B[i]])

        log10p = sum(result)
        log10p_prime = np.sum(result_prime, axis=0)

        delta_h = sum(enthalpy)
        enthalpy_prime = np.sum(enthalpy_prime, axis=0)

        p_sat = (10**log10p)*101325 # Pa

        print('at T = ', T, 'K')
        print('P_sat =', p_sat, ' Pa' )
        print('Delta_H =', delta_h, ' kJ/mol')
        print('**************************************')
        if Antoine == True:
            # Fitting process to find Antonie coefficients for Log P = A - (B / (T + C))
            T_min_values = 273.15
            T_max_values = 390.15
            T_values = np.linspace(T_min_values, T_max_values, 1000)
            P_values = log10p_prime[0]/T_values + log10p_prime[1] + log10p_prime[2]*T_values + log10p_prime[3]*np.log(T_values)

            # Define the function to fit
            def fit_function(T, A, B, C):
                return A - (B / (T + C))

            # Fit the data
            popt, pcov = curve_fit(fit_function, T_values, P_values, maxfev=500000)

            # Extract the fitted parameters
            A_fit, B_fit, C_fit = popt
        ###########################################################################
            if yieldABC:
                yield A_fit, B_fit, C_fit
        # write output line by line for each compound
        if line in inorg_list:
            output_line = (
                            f"-, -, {line_intact}, {column2}, {column7}, {column8}, {column9}, {column10}, "\
                            f"{column11}, {column12}, {column13}, {column14}, {column15}, {column16}, " \
                            f"{column17}, {column18}, {column19}, {column20}, {column21}, {column22}, " \
                            f"{column23}, {column24}, {column25}, {column26}, {column27}, {column28}, " \
                            f"{column29}, {column30}, {column31}, {column32}, {column33}, {column34}, " \
                            f"{column35}, " \
                            f"-, " \
                            f"-, -\n"
                            )

        else:
            if Antoine == True:
                output_line = (
                                f"{p_sat}, {delta_h}, {line_intact}, {column2}, {column7}, {column8}, {column9}, {column10}, "\
                                f"{column11}, {column12}, {column13}, {column14}, {column15}, {column16}, " \
                                f"{column17}, {column18}, {column19}, {column20}, {column21}, {column22}, " \
                                f"{column23}, {column24}, {column25}, {column26}, {column27}, {column28}, " \
                                f"{column29}, {column30}, {column31}, {column32}, {column33}, {column34}, " \
                                f"{column35}, " \
                                f"{log10p_prime[0]}/T + ({log10p_prime[1]}) + ({log10p_prime[2]}*T) + ({log10p_prime[3]}*log(T)), " \
                                f"{A_fit} - ({B_fit})/(T + {C_fit}), {enthalpy_prime[0]} - ({enthalpy_prime[2]}/T**2) - ({enthalpy_prime[3]}/T)\n"
                                )

            else:
                output_line = (
                                f"{p_sat}, {delta_h}, {line_intact}, {column2}, {column7}, {column8}, {column9}, {column10}, "\
                                f"{column11}, {column12}, {column13}, {column14}, {column15}, {column16}, " \
                                f"{column17}, {column18}, {column19}, {column20}, {column21}, {column22}, " \
                                f"{column23}, {column24}, {column25}, {column26}, {column27}, {column28}, " \
                                f"{column29}, {column30}, {column31}, {column32}, {column33}, {column34}, " \
                                f"{column35}, " \
                                f"{log10p_prime[0]}/T + ({log10p_prime[1]}) + ({log10p_prime[2]}*T) + ({log10p_prime[3]}*log(T)), " \
                                f"-, {enthalpy_prime[0]} - ({enthalpy_prime[2]}/T**2) - ({enthalpy_prime[3]}/T)\n"
                                )


        file_out.write(output_line)
        s_number = s_number + 1

    file_out.close()
    if isinstance(filename_in, str):
        file_in.close()

if __name__ == '__main__':
    filename_in = "SMILES.txt"  #  input filename
    filename_out = "output.txt"  # output filename
    simpol_saturation_pressure(filename_in, filename_out)

    # Read data from the text file
    with open('output.txt', 'r') as file:
        lines = file.readlines()

    # Split each line by commas and organize into rows and columns
    data = [line.strip().split(',') for line in lines]

    # Write the data to an Excel file
    with open('output.csv', 'w') as file:
        # Iterate over each line in the data
        for line in data:
            # Join the values with commas and write to the file
            file.write(','.join(line) + '\n')

    print('Done!!!')


################################################################
