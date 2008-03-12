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
_default_timescale = os.path.join(data_path, 'data', 'candekent.dat')


def get_assymetry(assymetry_file=None):
    """
    retrieves assymetry percentage and corresponding time.
    If no file is provided assymetry is assumed to be zero
    throughout. Returns a dictionary consisting of
    {start_of_period:assymetry_for_period}
    i.e. beginning of time period maps to assymetry of
    that particular period.
    """

    if assymetry_file is None: return {0:0}

    assymetry = {}
    
    file_content = open(assymetry_file).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of file:
        period_start period_end assymetry

        where
        period_start is start of period in myrs (million years)
        period_end   is end of period in myrs
        assymetry    is the assymetry for that period in %

        % at start of line is a comment
        
        Observe that period_end is only for the user's convenience
        and the only period_start is used to indicate the assymetry
        of a period until next period starts
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue

        columns = line.split()
        #columns[1] is start, columns[2] is assymetry
        assymetry[-eval(columns[1])] = {'assymetry':eval(columns[2])}

    return assymetry


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

    file_content = open(jump_file).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of file:
        time distance

        where
        time     is time of jump myrs (million years)
        distance is the distance of that jump in km

        % at start of line is a comment
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue

        columns = line.split()
        #columns[0] is time, columns[1] is distance
        jumps[-eval(columns[0])] = {'jump':eval(columns[1])}

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
    
    file_content = open(spreadingrate_file).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of file:
        period_start period_end spreading_rate

        where
        period_start   is start of period in myrs (million years)
        period_end     is end of period in myrs
        spreading_rate is the spreading rate for that period in km/myr

        % at start of line is a comment
        
        Observe that period_end is only for the user's convenience
        and the only period_start is used to indicate the spreading
        rate of a period until next period starts
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue
        
        columns = line.split()
        #columns[1] is start, columns[2] is full spreading rate (half needed)
        spr_rates[-eval(columns[1])] = {'spreadingrate':(eval(columns[2])/2)}

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
    
    file_content = open(timescale).read()
    lines = file_content.splitlines()

    number_of_periods = len(lines)
    if ((number_of_periods % 2) == 0): polarity = 'r'
    else: polarity = 'n'
    
    for line in lines:
        """
        Format of file:
        period_start period_end period_name

        where
        period_start is start of period in myrs (million years)
        period_end   is end of period in myrs
        period_name  is name of period

        % at start of line is a comment
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue
        
        columns = line.split()
        
        reversed_timescale[-eval(columns[1])] = {'polarity':polarity}

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

    file_content = open(input_file).read()
    lines = file_content.splitlines()
    for line in lines:
        """
        Format of file:
        distance longditude latitude depth anomaly

        where
        distance   is the distance in kilometers
        longditude is longditude coordinate (not needed)
        latidtude  is latitude coordinate (not needed)
        depth      is depth in kilometers from ocean top (must be negated)
        anomaly    is anomaly of magnetic measurements in nanoTesla

        % at start of line is a comment
        """

        #Ignore comments and blank lines
        if re.match('^(%)|(\s*$)',line):
            continue
        
        columns = line.split()
        distance.append(eval(columns[0])) #Evaluate as number
        depth.append(-1*eval(columns[3])) #Evaluate as number and negate
        anomaly.append(eval(columns[4])) #Evaluate as number

    return (distance, depth, anomaly)

