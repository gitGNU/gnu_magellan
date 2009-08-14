# -*- coding: utf-8 -*-

"""
data.py - gathers data from different files and returns data structures

Copyright (C) 2008 Tryggvi Björgvinsson <tryggvib@hi.is>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os,sys,re
import Magellan

# Filenames for default files
data_path = os.path.split(Magellan.__file__)[0]
_default_timescale = os.path.join(data_path, 'data', 'Lourens2004.dat')

def _read_period_file(period_file, key_column, value_column):
    """
    A generator function which reads information from files
    with time periods and a specific value for those periods
    Theoretically it can be used for similar, non-time period
    files since it returns the content of given columns (through
    the parameters key_column and value_column). Returns a tuple
    containing the evaluated value of the key column and the
    value column
    """

    filepath = os.path.expanduser(period_file)
    file_content = open(filepath).read()
    lines = file_content.splitlines()
    
    for line in lines:
        """
        Format of files:
        period_start period_end value

        where
        period_start is start of period in myrs (million years)
        period_end   is end of period in myrs
        value        is the value for that period

        % at start of line is a comment
        
        Observe that period_end is only for the user's convenience
        and the only period_start is used to indicate the asymmetry
        of a period until next period starts
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue

        #Split line into columns
        columns = line.split()

        yield (eval(columns[key_column]), eval(columns[value_column]))

def get_asymmetry(asymmetry_file=None, key_name='asymmetry'):
    """
    retrieves asymmetry percentage and corresponding time.
    Returns a dictionary consisting of
    {start_of_period:asymmetry_for_period}
    i.e. beginning of time period maps to asymmetry of
    that particular period.
    If no file is provided an empty directory is returned.
    """

    if asymmetry_file is None: return {}

    asymmetry = {}

    for (key,value) in _read_period_file(asymmetry_file,1,2):
        asymmetry[-key] = {key_name:value}

    return asymmetry

def get_magnetization(magnetization_file=None):
    """
    retrieves magnitude of magnetization and corresponding time.
    Returns a dictionary consisting of
    {start_of_period:magnetization_for_period}
    i.e. beginning of time period maps to magnetization of
    that particular period.
    If no file is provided an empty directory is returned
    """

    if magnetization_file is None: return {}

    # It is actually the same process as with asymmetry, so we use
    # the asymmetry function
    return get_asymmetry(magnetization_file,'magnetization')

def get_jumps(jump_file=None):
    """
    retrieves jump distances and corresponding time.
    If no file is provided it is assumed that no jumps
    have taken place. Returns a dictionary consisting of
    {time_of_jump:jump_distance}
    i.e. time of jump maps to distance of that particular
    jump.
    """

    if jump_file is None: return {}

    jumps = {}

    for (key,value) in _read_period_file(jump_file,0,1):
        jumps[-key] = {'jump':value}
   
    return jumps


def get_spreadingrate(spreadingrate_file):
    """
    retrieves spreading rates and corresponding time. A
    filename must be provided. Returns a dictionary
    consisting of {start_of_period:spreading_rate_for_period},
    i.e. beginning of time period maps to spreading rate of
    that particular period.
    """

    spr_rates = {}
    
    for (key,value) in _read_period_file(spreadingrate_file,1,2):
        spr_rates[-key] = {'spreadingrate':(value/2)}

    return spr_rates


def get_timescale(timescale=None):
    """
    gets the reversed timescale. If no input file is
    provided, use default file of Cande Kent from '95.
    Assumes normal polarity at time zero.
    Returns reversed timescale as a list of tuples:
    [(name of period, start of period, end of period)...]
    """
    
    if timescale is None: timescale = _default_timescale

    reversed_timescale = {}
    
    filepath = os.path.expanduser(timescale)
    file_content = open(filepath).read()
    lines = file_content.splitlines()

    number_of_periods = len(lines)
    #if ((number_of_periods % 2) == 0): polarity = 'r'
    #else: polarity = 'n'
    polarity='n'
    for (key,value) in _read_period_file(filepath,1,0):
        reversed_timescale[-key] = {'polarity':polarity}

        #Swap polarities
        if polarity == 'n': polarity = 'r'
        else: polarity = 'n'
    
    return reversed_timescale


def get_trackdata(input_file):
    """
    gathers data from track file. Input file must be
    provided. Data gathered is distance, anomaly and
    depth. Returns...
    """
    
    distance = []
    anomaly = []
    depth = []
    distance = []
    depth = []

    filepath = os.path.expanduser(input_file)
    file_content = open(filepath).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of file:
        distance longtitude latitude depth anomaly

        where
        distance   is the distance in kilometers
        longtitude is longditude coordinate (not needed)
        latitude  is latitude coordinate (not needed)
        depth      is depth in kilometers from ocean top (must be negated)
        anomaly    is anomaly of magnetic measurements in nanoTesla

        % at start of line is a comment
        """

        # Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue
        
        columns = line.split()
        distance.append(eval(columns[0])) # Evaluate as number
        depth.append(eval(columns[3])) # Evaluate as number and negate
        anomaly.append(eval(columns[4])) # Evaluate as number

    if distance[0] > distance[1]:
	distance.reverse()
	depth.reverse()
	anomaly.reverse()
    # Extending the bathymetry 20 kms in both directions to avoid edge affects.
    depth_calc = [depth[0]]*20 + depth + [depth[-1]]*20
    distance_first = [distance[0] - x for x in range(1,21)]
    distance_first.reverse()
    distance_calc = distance_first + distance + [x + distance[-1] for x in range(1,21)]
   
    return (distance_calc, depth_calc, distance, anomaly)

def get_configurations(config_file=None):
    """
    Go through a configuration file (project file)
    where the user can set different values to settings
    otherwised provided through flags. Makes it easier
    to manage flags and use magellan
    Returns an empty dictionary if no file is provided
    """

    configurations = {}

    if config_file is None: return configurations

    file_content = open(config_file).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of files:
        item=value

        where
        item is the configurable setting
        value is the value for that item

        % at start of line is a comment
        
        Configurations which do something in magellan are:
        
        asymmetry = location of asymmetry file
        data = location of data file (track data)
        jump = location of jumps file
        magnetization = location of magnetization file
        spreading rate = location of spreading rate file
        
        declination = amount of declination
        inclination = amount of inclination
        obliquity = angle of obliquity
	azimuth = azimuth of the ridge relative to north
        thickness = thickness of magnetized layer

        graphs = which graphs to plot (not implemented yet)
        """

        # Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue

        # Split line into columns by '='
        # Remove preceeding and trailing whitespaces
        # Pick out the words on each side of 
        pattern = re.compile(r'^\s*(\w.*\w)\s*=\s*(\S.*\S)\s*$')
        match = pattern.match(line)
        
        configurations[match.group(1).lower()] = match.group(2)

    return configurations
