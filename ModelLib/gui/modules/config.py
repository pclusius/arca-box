#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=============================================================================
Copyright (C) 2022 Multi-Scale Modelling group
Institute for Atmospheric and Earth System Research (INAR),
University of Helsinki
Contact information sosaa@helsinki.fi

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

import configparser
import pathlib

from .resources import resource_path, open_resource

_CONFIG_PATH = "conf/config.ini"


# Get a configuration value from the INI file
def get_config(*args, **kwargs):
    gui_config = configparser.ConfigParser()

    if pathlib.Path(resource_path(_CONFIG_PATH)).exists():
        with open_resource(_CONFIG_PATH, "r") as file:
            gui_config.read_file(file)

    return gui_config.get(*args, **kwargs)


# Write a configuration value to the INI file
def set_config(section, option, value):
    gui_config = configparser.ConfigParser()

    if pathlib.Path(resource_path(_CONFIG_PATH)).exists():
        with open_resource(_CONFIG_PATH, "r") as file:
            gui_config.read_file(file)

    if not gui_config.has_section(section):
        gui_config.add_section(section)

    gui_config.set(section, option, value)

    with open_resource(_CONFIG_PATH, "w") as file:
        gui_config.write(file)


# Remove a configuration option from section in the INI file
def remove_config(section, option):
    if not pathlib.Path(resource_path(_CONFIG_PATH)).exists():
        return

    gui_config = configparser.ConfigParser()

    with open_resource(_CONFIG_PATH, "r") as file:
        gui_config.read_file(file)

    if gui_config.has_option(section, option):
        gui_config.remove_option(section, option)

    with open_resource(_CONFIG_PATH, "w") as file:
        gui_config.write(file)
