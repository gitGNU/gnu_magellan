# -*- coding: utf-8 -*-

"""
calc.py - core computations for magellan.

Copyright (C) 2008  Tryggvi Björgvinsson <tryggvib@hi.is>

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

from math import cos, sin, atan2, radians, degrees, sqrt, log

_default_thickness = '0.5'
_default_inclination = '75'
_default_declination = '18'

def create_change_timeline(assym,spread,jump,time):
    """
    creates a timeline of changes from arrays containing
    assymetry, spreading rates, jumps, and timescale for
    easier lookup when processing bathymetry. Returns a
    list of tuples (each tuple contains start of period
    and dictionary of changing values) sorted by start
    of period:
    [(start_of_period, {change:value})]
    """

    assym_spread = set(assym).intersection(set(spread))
    for same in assym_spread:
        spread[same].update(assym[same])

    assym.update(spread)

    assym_jump = set(assym).intersection(set(jump))
    for same in assym_jump:
        jump[same].update(assym[same])

    assym.update(jump)

    assym_time = set(assym).intersection(set(time))
    for same in assym_time:
        time[same].update(assym[same])

    assym.update(time)

    assym[0] = {'assymetry':0, 'spreadingrate':0}

    return [(k,assym[k]) for k in sorted(assym.keys())]

# 12.july 2008
# Here we need to add the polarity
def create_anomaly_model(dist,deep,parameters):
    """
    creates an anomaly model from distance and depth.
    Other parameters needed are thickness, declination,
    inclination, azimuth, and magnetizsation
    The model is based on theoretical computations.
    Returns a list of depths sorted by distance in x
    direction: [depth]
    """

    thickness = eval(parameters.get('thickness', _default_thickness))
    declination = radians(eval(parameters.pop('declination',
                                              _default_declination)))
    inclination = radians(eval(parameters.pop('inclination',
                                              _default_inclination)))

    dmy_azimuth = radians(100) # Will be given by parameter file
    dmy_magnetization = 150 # Will be given by parameter file

    sinI = sin(inclination)
    cosI = cos(inclination)
    cosCminD = cos(dmy_azimuth-declination)

    Jx = dmy_magnetization * cosI * cosCminD
    Jz = dmy_magnetization * sinI

    model = {}
    for distance in dist:
        model[distance] = 0

    z1 = deep[0]
    x1 = dist[0]
    for position in (range(len(dist))[1:]):

        z2 = deep[position]
        x2 = dist[position]

        p_former = ((z2-z1)**2)/((z2-z1)**2 + (x1-x2)**2)
        p_latter = ((z2-z1)*(x1-x2))/((z2-z1)**2 + (x1-x2)**2)

        for distance in dist:
            x1_d = (x1-distance)**2
            x2_d = (x2-distance)**2

            theta1 = degrees(atan2(z1, x1_d))
            theta2 = degrees(atan2(z2, x2_d))

            r2_1 = sqrt((x2_d + z2**2)) / sqrt((x1_d + z1**2))

            P = p_former * (theta1-theta2) + p_latter * 0.5 * log(r2_1)
            Q = p_latter * (theta1-theta2) - p_former * 0.5 * log(r2_1)

            V = 2*(Jx*Q - Jz*P)
            H = 2*(Jx*P + Jz*Q)
            T = V*sinI + \
                H*cosI*cosCminD
            
            model[distance] += T

            theta1 = degrees(atan2((z1-thickness), x1_d))
            theta2 = degrees(atan2((z2-thickness), x2_d))

            r2_1 = (sqrt(x2_d + (z2-thickness)**2)) / \
                   (sqrt(x1_d + (z1-thickness)**2))

            P = p_former * (theta1-theta2) + p_latter * 0.5 * log(r2_1)
            Q = p_latter * (theta1-theta2) - p_former * 0.5 * log(r2_1)

            V = 2*(Jx*Q - Jz*P)
            H = 2*(Jx*P + Jz*Q)

            T = V*sinI + \
                H*cosI*cosCminD
            
            model[distance] += T

        z1 = z2
        x1 = x2

    return [model[k] for k in sorted(model.keys())]


def create_deltax(timeline):
    """
    Skil þetta ekki alveg
    'create the total difference traveled between
    changes in a timeline given a change timeline
    as input'. Returns a tuple of lists containing
    a tuple of distance and a dictionary which
    contains polarity, pseudo-faults and failed
    rifts:
    ([(distance, {polarity,pseudo faults,failed rift})],
     [(distance, {polarity,pseudo faults,failed rift})])

    The former tuple is distances in left direction and
    the latter tuple is distances in right direction

    polarity is either the string 'n' (normal) or 'r' (reverse)
    pseudo_fault and failed rift are booleans True if it is either
    """

    # Delta movement in right direction
    deltax_r = []
    # Delta movement in left direction
    deltax_l = []

    # Initialize with time zero
    (prev_time, action) = timeline.pop(0)

    assymetry = 0 # Default assymetry
    if action.has_key('assymetry'):
        assymetry = action['assymetry']
    
    spreading_rate = action['spreadingrate']
    polarity = action['polarity']

    pseudo_fault = False # Let know of pseudo faults
    
    # Spreading rate in right direction
    spread_r = spreading_rate * (1 + assymetry)
    # Spreading rate in left direction
    spread_l = spreading_rate * (1 - assymetry)

    for (time,action) in timeline:

        delta_t = time - prev_time
        # Delta distance in right direction
        distance_r = delta_t * spread_r
        # Delta distance in left direction
        distance_l = -delta_t * spread_l
        
        deltax_r.append((distance_r, {'polarity':polarity,
                                      'pseudo fault':pseudo_fault,
                                      'failed rift':False}))
        deltax_l.append((distance_l, {'polarity':polarity,
                                      'pseudo fault':pseudo_fault,
                                      'failed rift':False}))

        pseudo_fault = False

        if action.has_key('polarity'):
            polarity = action['polarity']

        if action.has_key('spreadingrate'):
            spreading_rate = action['spreadingrate']
            spread_r = spreading_rate * (1 + assymetry)
            spread_l = spreading_rate * (1 - assymetry)

        if action.has_key('assymetry'):
            assymetry = action['assymetry']
            spread_r = spreading_rate * (1 + assymetry)
            spread_l = spreading_rate * (1 - assymetry)

        # Jump computations
        if action.has_key('jump'):
            jump = action['jump']
            if (jump > 0):
                (deltax_r, deltax_l) = _jump_reposition(deltax_r,deltax_l,jump)
            if (jump < 0):
                (deltax_l, deltax_r) = _jump_reposition(deltax_l,deltax_r,jump)

            pseudo_fault = True

        prev_time = time

    deltax_l.reverse()
    deltax_r.reverse()
    return (deltax_l, deltax_r)

def _jump_reposition(jump_in, move_to, jump):
    """
    repositions the ridge axis when a jump occurs.
    The method takes in the jump_in which represents
    the half which the jump will land in and moves
    every distance up to the jump point into the
    move_to half, setting pseudo faults and failed rifts
    accordingly. Returns the recomputed halfs in a
    tuple: (jump_in, move_to)
    """

    pseudo_fault = False
    failed_rift = True
    (dx,prefs) = jump_in.pop()

    while(jump > dx):
        fr_tmp = failed_rift
        failed_rift = prefs['failed rift']
        prefs['failed rift'] = fr_tmp

        pf_tmp = pseudo_fault
        pseudo_fault = prefs['pseudo fault']
        prefs['pseudo fault'] = pf_tmp

        move_to.append((-dx,prefs))
        jump -= dx

        (dx,prefs) = jump_in.pop()

    new_prefs = {}
    new_prefs['polarity'] = prefs['polarity']
    new_prefs['failed rift'] = failed_rift
    new_prefs['pseudo fault'] = pseudo_fault

    move_to.append((-jump, new_prefs))
    jump_in.append(((dx-jump),prefs))

    return(jump_in,move_to)
        

def create_magnetized_layer(deltax_l, deltax_r, min_l, max_r):
    """
    create a magnetized layer between min_l (maximum distance
    to the left - or minimum distance since it is a negative
    number) and max_r (maximum distance to the right) from the
    two halfs. Returns a list of tuples. Each tuple includes
    another tuple with start and end postitions, and the
    polarity for that distance:
    [((start,end),polarity)]
    """

    magnetized_layer = []
    distance_sum = 0
    start = 0
    polarity = 'n'

    sum_right = []
    sum_left = []

    for (deltax, point_prefs) in deltax_r:

        if (distance_sum > max_r):
            sum_right.append(((start, max_r), polarity))
            break

        if (polarity != point_prefs['polarity']):
            sum_right.append(((start, round(distance_sum)), polarity))

            start = round(distance_sum)
            polarity = point_prefs['polarity']
            
        distance_sum += deltax

    distance_sum = 0
    start = 0
    polarity = 'n'
    
    for (deltax, point_prefs) in deltax_l:

        if(distance_sum < min_l):
            sum_left.append(((min_l,start), polarity))
            break

        if (polarity != point_prefs['polarity']):
            sum_left.append(((round(distance_sum), start), polarity))

            start = round(distance_sum)
            polarity = point_prefs['polarity']
            
        distance_sum += deltax

    sum_left.reverse()

    ((start_l,end_l),polarity) = sum_left[-1]
    ((start_r,end_r),polarity) = sum_right[0]

    magnetized_layer.extend(sum_left[:-1])
    magnetized_layer.append(((start_l,end_r),polarity))
    magnetized_layer.extend(sum_right[1:])

    return magnetized_layer

def create_faults_and_rifts(deltax_l, deltax_r, min_l, max_r):
    """
    
    """

    faults_and_rifts = []

    distance_sum = 0

    for (deltax, point_prefs) in deltax_l:

        if (distance_sum < min_l): break

        if (point_prefs['pseudo fault'] or point_prefs['failed rift']):
            faults_and_rifts.append((round(distance_sum),
                                     point_prefs['pseudo fault'],
                                     point_prefs['failed rift']))
    
        distance_sum += deltax


    distance_sum = 0

    for (deltax, point_prefs) in deltax_r:

        if (distance_sum > max_r): break

        if (point_prefs['pseudo fault'] or point_prefs['failed rift']):
            faults_and_rifts.append((round(distance_sum),
                                     point_prefs['pseudo fault'],
                                     point_prefs['failed rift']))
            
        distance_sum += deltax

    return faults_and_rifts
