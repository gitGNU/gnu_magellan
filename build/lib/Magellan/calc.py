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

from math import cos, sin, atan2, radians, degrees, sqrt, log, pi

_default_thickness = '0.5'
 # This has to be a decimal number
_default_azimuth = '90'
_default_inclination =  '45'
_default_declination = '45'
_default_obliquity = '0'
# The obliquity variable has to be global, it is used in two defs
obliquity = 0

def create_change_timeline(assym,spread,jump,magnet,time):
    """
    creates a timeline of changes from arrays containing
    assymetry, spreading rates, jumps, and timescale for
    easier lookup when processing bathymetry. Returns a
    list of tuples (each tuple contains start of period
    and dictionary of changing values) sorted by start
    of period:
    [(start_of_period, {change:value})]
    Note: Updates the asymmetry dictionary
    """
    
    assym_spread = set(assym).intersection(set(spread))
    for same in assym_spread:
        spread[same].update(assym[same])

    assym.update(spread)

    assym_jump = set(assym).intersection(set(jump))
    for same in assym_jump:
        jump[same].update(assym[same])

    assym.update(jump)

    assym_magnet = set(assym).intersection(set(magnet))
    for same in assym_magnet:
        magnet[same].update(assym[same])

    assym.update(magnet)

    assym_time = set(assym).intersection(set(time))
    for same in assym_time:
        time[same].update(assym[same])

    assym.update(time)

    assym[0] = {'assymetry':0, 'spreadingrate':0}

    return [(k,assym[k]) for k in sorted(assym.keys())]

def create_anomaly_model(dist,deep,parameters, magnet_layer):
    """
    creates an anomaly model from distance and depth.
    Other parameters needed are thickness, declination,
    inclination, azimuth, and magnetizsation
    The model is based on theoretical computations.
    Returns a list of depths sorted by distance in x
    direction: [depth]
    """

    thickness = eval(parameters.pop('thickness', _default_thickness))
    declination = radians(eval(parameters.pop('declination', _default_declination)))
   
    inclination = radians(eval(parameters.pop('inclination', _default_inclination)))
    global obliquity  
    
    azimuth = radians(eval(parameters.pop('azimuth', _default_azimuth)))
    print obliquity
   
    print thickness, inclination, declination, obliquity, azimuth
    dmy_magnetization = 2 # Will be given by parameter file
    
    # Here we have to multiply with 4pi because we are working in the SI system but these equations were 'derived' 
    # for the cgs system. Basically k_cgs = 4pi k_si
    
    #magnetization_1 = 10 # * 4 * pi # = k * H_e = susceptibility * scalar_earth_magnetic_field_strenght = magnetization
    sus = 0.001
    #magnetization = magnetization_1/(sus*4*pi)
    mag_dict = {}
    projected_dist = []
    for distance in dist:
    	projected_dist.append(distance*cos(obliquity))

    for ((first_pos, last_pos),pol,magnet) in magnet_layer:
	mag_dict[projected_dist[next_index(projected_dist,first_pos)]] = (pol,magnet)	

    #print mag_dict
    sinI = sin(inclination)
    cosI = cos(inclination)
    cosCminD = cos(azimuth-declination)
	
    pol_direction = 1


    model = {}
    for distance in projected_dist:
        model[distance] = 0
    
    # 1. Method taking all points in the magnetic layer
    
    z1 = deep[0] 
    x1 = projected_dist[0]
    for position in range(1,len(projected_dist)):
	z2 = deep[position]
	x2 = projected_dist[position]
	if x1 in mag_dict:
	    if (mag_dict[x1][0] == 'n'):
	        pol_direction = 1
	    else:
	        pol_direction = -1
	    mag_field =  mag_dict[x1][1]/(sus*4*pi)
    

	z3 = z1 + thickness
	z4 = z2 + thickness
	z1_pow2 = z1**2    
	z2_pow2 = z2**2    
	z3_pow2 = z3**2
	z4_pow2 = z4**2
       
	"""
    # 2.Method calculating the contribution from the magnetic blocks only 
    for ((x1, x2),pol) in magnet_layer:
	if x1 == x2:
	    continue

	position_x1 = next_index(dist,x1)
	z1 = deep[position_x1]
  
	if (pol == 'n'):
	    pol_direction = 1
	else:
	    pol_direction = -1


	position_x2 = next_index(dist, x2)
        z2 = deep	[position_x2]
	z3 = z1 + thickness
	z4 = z2 + thickness
	z1_pow2 = z1**2    
	z2_pow2 = z2**2    
	z3_pow2 = z3**2
	z4_pow2 = z4**2
        """  
        
	#### Talwani or Won & Bevis ####    
	# Talwani
	Jx = mag_field*cosI*cosCminD
    	Jz = mag_field*sinI
        for distance in projected_dist:
            x1_calc = x1 - distance
    	    x2_calc = x2 - distance

	    theta1 = atan2(z2,x2_calc)
	    theta2 = atan2(z4,x2_calc)
	    theta3 = atan2(z3,x1_calc)
	    theta4 = atan2(z1,x1_calc)

	    x2_calc_pow2 = x2_calc**2
	    x1_calc_pow2 = x1_calc**2
	    # Right surface; from (x2,z2) to (x2,z4)
	    z21 = z4-z2
	    x12 = 0
            r1 = sqrt(x2_calc_pow2 + z2_pow2)
	    r2 = sqrt(x2_calc_pow2 + z4_pow2)
	    #print "r1:",r1,"r2:",r2, "x1:", x1, "x2_calc_pow:",x2_calc_pow2
	    P_r = (theta1-theta2)
	    Q_r = -1*log(r2/r1)

	    V_r = 2*(Jx*Q_r - Jz*P_r)
	    H_r = 2*(Jx*P_r + Jz*Q_r)
	    T_r = V_r*sinI + H_r*cosI*cosCminD

	    # Left surface; from (x1,z3) to (x1,z1)
	    z21 = z1-z3
	    x12 = 0
	    r1 = sqrt(x1_calc_pow2 + z4_pow2)
	    r2 = sqrt(x1_calc_pow2 + z1_pow2)

	    P_l = (theta3-theta4)
	    Q_l = -1*log(r2/r1)

	    V_l = 2*(Jx*Q_l - Jz*P_l)
	    H_l = 2*(Jx*P_l + Jz*Q_l)
	    T_l = V_l*sinI + H_l*cosI*cosCminD
	    
	    # Top surface; from (x1,z1) to (x2,z2)
	    z21 = z2-z1
	    x12 = x1 - x2
	    r1 = sqrt(x1_calc_pow2 + z1_pow2)
	    r2 = sqrt(x2_calc_pow2 + z2_pow2)
	
	    const1 = z21**2/(z21**2 + x12**2)
	    const2 = z21*x12/(z21**2 + x12**2)
	    P_t = const1*(theta4 - theta1) + const2*log(r1/r2)
	    Q_t = const2*(theta4-theta1) - const1*log(r1/r2)

	    V_t = 2*(Jx*Q_t - Jz*P_t)
	    H_t = 2*(Jx*P_t + Jz*Q_t)
	    T_t = V_t*sinI + H_t*cosI*cosCminD

	    # Bottom surface; from (x2,z4) to (x1,z3)
	    z21 = z3-z4
	    x12 = x2-x1
	    r1 = sqrt(x2_calc_pow2 + z4_pow2)
	    r2 = sqrt(x1_calc_pow2 + z3_pow2)
	
	    const1 = z21**2/(z21**2 + x12**2)
	    const2 = z21*x12/(z21**2 + x12**2)
	    P_b = const1*(theta2-theta3) + const2*log(r1/r2)
	    Q_b = const2*(theta2-theta3) - const1*log(r1/r2)

	    V_b = 2*(Jx*Q_b - Jz*P_b)
	    H_b = 2*(Jx*P_b + Jz*Q_b)
	    T_b = V_b*sinI + H_b*cosI*cosCminD
            model[(distance)] += pol_direction * (T_b + T_t + T_l + T_r)

	z1 = z2
	x1 = x2 
    return [model[k] for k in sorted(model.keys())]
    """

	# Won and Bevis
        R1 =  (x2 - x1)**2 + (z2 - z1)**2
        R2 =  (z4 - z2)**2 # (x2 - x2)**2 + (z4 - z2)**2
        R3 =  (x2 - x1)**2 + (z3 - z4)**2
        R4 =  (z3 - z1)**2 # (x1 - x1)**2 + (z3 - z1)**2


        z21 = z2 - z1
        z42 = z4 - z2
        z34 = z3 - z4
        z13 = z1 - z3
        x21 = x2 - x1
        x12 = x1 - x2

	# NB constant2 = constant4 = 0
	constant1 = (x2 - x1) / R1
	constant3 = (x2 - x1) / R3

        
        for distance in dist:
	    distance_calc = distance
            x1_d = (x1-distance_calc)**2
            x2_d = (x2-distance_calc)**2
	    xd1 = x1 - distance_calc
	    xd2 = x2 - distance_calc
	    # Added + for z1 - thickness. esThis is a right handed coordinate system like Talwani describes. +Z axis 		      points down and the +X axis east and +Y axis south.
	    theta1 = atan2(z1, xd1)
            theta2 = atan2(z2, xd2)

            theta3 = atan2(z3, xd1)
            theta4 = atan2(z4, xd2)

            r1_21 = 1.0*(x2_d + z2**2) / (x1_d + z1**2)
            r2_21 = 1.0*(x2_d + z4**2) / (x2_d + z2**2)
            r3_21 = 1.0*(x1_d + z3**2) / (x2_d + z4**2)
            r4_21 = 1.0*(x1_d + z1**2) / (x1_d + z3**2)
	
	    P1 = ((xd1*z2-xd2*z1) / R1) * ((xd1*x21 - z1*z21) / (x1_d + z1**2) - (xd2*x21 - z2*z21) / (x2_d+z2**2))
	    Q1 = ((xd1*z2-xd2*z1) / R1) * ((xd1*z21 + z1*x21) / (x1_d + z1**2) - (xd2*z21 + z2*x21) / (x2_d+z2**2))

	    P2 = ((xd2*z4-xd2*z2) / R2) * (( -1 * z2 * z42) / (x2_d + z2**2) - (-1* z4*z42) / (x2_d + z4**2))
	    Q2 = ((xd2*z4-xd2*z2) / R2) * (xd2*z42 / (x2_d + z2**2) - xd2*z42 / (x2_d + z4**2))

	    P3 = ((xd2*z3-xd1*z4) / R3) * ((xd2*x12 - z4*z34) / (x2_d + z4**2) - (xd1*x12 - z3*z34) / (x1_d+z3**2))
	    Q3 = ((xd2*z3-xd1*z4) / R3) * ((xd2*z34 + z4*x12) / (x2_d + z4**2) - (xd1*z34 + z3*x12) / (x1_d+z3**2))

	    P4 = ((xd1*z1-xd1*z3) / R4) * (( -1 * z3*z13) / (x1_d + z3**2) - (-1 * z1*z13) / (x1_d + z1**2))
	    Q4 = ((xd1*z1-xd1*z3) / R4) * (xd1*z13 / (x1_d + z3**2) - xd1*z13 / (x1_d + z1**2))

	   
	    theta12 = theta1 - theta2
	    theta24 = theta2 - theta4
	    theta43 = theta4 - theta3
	    theta31 = theta3 - theta1
	    
	    dZdz1 = constant1 * (x21 * theta12 + z21 * 0.5 * log(r1_21)) - P1
	    dZdx1 = -1 * constant1 * (z21 * theta12 + z21**2 * 0.5 * log(r1_21)) + Q1
	    dXdz1 = -1 * constant1 * (z21 * theta12 - x21 * 0.5 * log(r1_21)) + Q1
	    dXdx1 = constant1 * ((z21**2 / x21) * theta12 - z21 * 0.5 * log(r1_21)) + P1

	    dZdz2 = -1 *  P2
	    dZdx2 = Q2
	    #dZdx2 = -1 * z42**2 / R2 * 0.5 * log(r2_21) + Q2
	    dXdz2 = Q2
	    dXdx2 = P2
	    #dXdx2 = z42**2 / R2 * theta24 + P2

	    dZdz3 = constant3 * (x12 * theta43 + z34 * 0.5 * log(r3_21)) - P3
	    dZdx3 = -1 * constant3 * (z34 * theta43 + z34**2 * 0.5 * log(r3_21)) + Q3
	    dXdz3 = -1 * constant3 * (z34 * theta43 - x12 * 0.5 * log(r3_21)) + Q3
	    dXdx3 = constant3 * ((z34**2 / x12) * theta43 - z34 * 0.5 * log(r3_21)) + P3

	    dZdz4 = -1 * P4
	    dZdx4 = Q4
	    #dZdx4 =  -1 * z13**2 / R4 * 0.5 * log(r4_21) + Q4
	    dXdz4 = Q4
	    dXdx4 = P4
	    #dXdx4 = z13**2 / R4 * theta31 + P4

            Hz1 = 2*mag_field*(sinI * dZdz1 + cosI * cosCminD * dZdx1)
	    Hx1 = 2*mag_field*(sinI * dXdz1 + cosI * cosCminD * dXdx1)
            H1 = Hz1*sinI +  Hx1*cosI*cosCminD
            
            Hz2 = 2*mag_field*(sinI * dZdz2 + cosI * cosCminD * dZdx2)
            Hx2 = 2*mag_field*(sinI * dXdz2 + cosI * cosCminD * dXdx2)
            H2 = Hz2*sinI + Hx2*cosI*cosCminD
           
            Hz3 = 2*mag_field*(sinI * dZdz3 + cosI * cosCminD * dZdx3)
            Hx3 = 2*mag_field*(sinI * dXdz3 + cosI * cosCminD * dXdx3)
            H3 = Hz3*sinI + Hx3*cosI*cosCminD

            Hz4 = 2*mag_field*(sinI * dZdz4 + cosI * cosCminD * dZdx4)
            Hx4 = 2*mag_field*(sinI * dXdz4 + cosI * cosCminD * dXdx4)
            H4 = Hz4*sinI + Hx4*cosI*cosCminD
            model[(distance)] += pol_direction * (H1 + H2 + H3 + H4)
	     
        z1 = z2
        x1 = x2 
   
    return [model[k] for k in sorted(model.keys())]
    """

def inv_project_anomaly_model(anomaly_model):
    """
    projects the anomaly_model back to the original track.
    """
    global obliquity
    print "hallo",obliquity
    projected_anomaly_model = []
    
    for value in anomaly_model:
    	projected_anomaly_model.append(value/cos(obliquity))

    return projected_anomaly_model

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
    
    magnetization = 10 #A/m -- # Default magnetization
    if action.has_key('magnetization'):
    	magnetization = action['magnetization']

    pseudo_fault = False # Let know of pseudo faults
    
    # Spreading rate in right direction
    spread_r = spreading_rate * (1 + assymetry)
    # Spreading rate in left direction
    spread_l = spreading_rate * (1 - assymetry)


    for (time,action) in timeline:

        delta_t = time - prev_time
        # Delta distance in right direction
        distance_r = delta_t * spread_r#/cos(radians(theta))
        # Delta distance in left direction
        distance_l = -delta_t * spread_l#/cos(radians(theta))
        
        deltax_r.append((distance_r, {'polarity':polarity,
        			      'magnetization':magnetization,
                                      'pseudo fault':pseudo_fault,
                                      'failed rift':False}))
        deltax_l.append((distance_l, {'polarity':polarity,
                		      'magnetization':magnetization,
                                      'pseudo fault':pseudo_fault,
                                      'failed rift':False}))

        pseudo_fault = False

        if action.has_key('polarity'):
            polarity = action['polarity']
            
        if action.has_key('magnetization'):
            magnetization = action['magnetization']

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

    while(abs(jump) > abs(dx)):
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
    new_prefs['magnetization'] = prefs['magnetization']
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
    magnet = 10 # Default magnetization

    sum_right = []
    sum_left = []

    for (deltax, point_prefs) in deltax_r:

        if (distance_sum > max_r):
            sum_right.append(((start, max_r), polarity, magnet))
            break

        if (polarity != point_prefs['polarity'] or magnet != point_prefs['magnetization']):
            sum_right.append(((start, distance_sum), polarity, magnet))

            start = distance_sum
            polarity = point_prefs['polarity']
            magnet = point_prefs['magnetization']
            
        distance_sum += deltax

    distance_sum = 0
    start = 0
    polarity = 'n'
    
    for (deltax, point_prefs) in deltax_l:

        if(distance_sum < min_l):
            sum_left.append(((min_l,start), polarity, magnet))
            break

        if (polarity != point_prefs['polarity'] or magnet != point_prefs['magnetization']):
            sum_left.append(((distance_sum, start), polarity, magnet))

            start = distance_sum
            polarity = point_prefs['polarity']
            magnet = point_prefs['magnetization']
            
        distance_sum += deltax

    sum_left.reverse()

    ((start_l,end_l),polarity,magnet) = sum_left[-1]
    ((start_r,end_r),polarity,magnet) = sum_right[0]

    magnetized_layer.extend(sum_left[:-1])
    magnetized_layer.append(((start_l,end_r),polarity,magnet))
    magnetized_layer.extend(sum_right[1:])

    return magnetized_layer

def create_projected_magnetized_layer(magnetized_layer,parameters):
    """
    projects previously made magnetized layer onto a profile 
    perpendicular to the ridge, if the original profile is obliuqe.
    Angle of projection is 90 - abs(ridge_orientation - track_orientation).
    """
    global obliquity
    obliquity = radians(eval(parameters.pop('obliquity', _default_obliquity)))
    print obliquity
    projected_magnetized_layer = []
    for ((start,end),polarity, magnet) in magnetized_layer:
	new_start = start * cos(obliquity)
	new_end = end*cos(obliquity)
	projected_magnetized_layer.append(((new_start,new_end),polarity,magnet))

    return projected_magnetized_layer


def create_faults_and_rifts(deltax_l, deltax_r, min_l, max_r):
    """
    
    """

    faults_and_rifts = []

    distance_sum = 0

    for (deltax, point_prefs) in deltax_l:

        if (distance_sum < min_l): break
	
	distance_sum += deltax
        if (point_prefs['pseudo fault'] or point_prefs['failed rift']):
            faults_and_rifts.append((distance_sum,
                                     point_prefs['pseudo fault'],
                                     point_prefs['failed rift']))

    distance_sum = 0

    for (deltax, point_prefs) in deltax_r:

        if (distance_sum > max_r): break

        distance_sum += deltax
        if (point_prefs['pseudo fault'] or point_prefs['failed rift']):
            faults_and_rifts.append((distance_sum,
                                     point_prefs['pseudo fault'],
                                     point_prefs['failed rift']))

    return faults_and_rifts

def next_index(indexed_list, search_value):
    for index in range(len(indexed_list)):
	if index == len(indexed_list)-1:
	    return index

	if (search_value < indexed_list[index+1]) and (search_value >= indexed_list[index]):
	    if (indexed_list[index+1] - search_value) > (search_value - indexed_list[index]):
		return index
	    else:
		return index+1

