# -*- coding: utf-8 -*-

"""
plot.py - plots data for magellan

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

from pylab import *
import Magellan

_default_thickness = '0.5'

def create_plot(dist, deep, anom, layer, faultrift, model, parameters):
    """
    Plot bathymetry profiles from distance, depth,
    anomalies, the magnetic layer and a model.
    Uses matplotlib to plot a nice graph.
    """

    thickness = eval(parameters.get('thickness', _default_thickness))
        
    fig = figure(figsize=(12,8))
    anomplot = fig.add_subplot(211)
    bathplot = fig.add_subplot(212)
    
    anomplot.set_title('Anomalies')
    anomplot.set_ylabel('nT')
    anomplot.plot(dist, anom)
    anomplot.plot(dist, model)
    anomplot.set_xlim(min(dist),max(dist))
    bathplot.set_title('Bathymetry')
    bathplot.set_xlabel('km')
    bathplot.set_ylabel('km')

    deepthick = map(lambda x: x-thickness, deep)

    for ((start,end),polarity) in layer:
        if polarity == 'n': fillcolor = 'b'
        else: fillcolor = 'w'
        
        index_lower = dist.index(start)
        index_upper = dist.index(end)
        
        if (index_lower != index_upper):
            xs = dist[index_lower:index_upper+1]
            lys = deep[index_lower:index_upper+1]
            tys = deepthick[index_lower:index_upper+1]

            x = concatenate( (xs, xs[::-1]) )
            y = concatenate( (lys, tys[::-1]) )
  
            bathplot.fill(x, y, facecolor=fillcolor)

    for (position,fault,rift) in faultrift:
        if fault:
            axvspan(position-0.5,position+0.5,facecolor='y')
        if rift:
            axvspan(position-0.5,position+0.5,facecolor='r')

    bathplot.plot(dist, deep, linewidth=2)
    bathplot.set_xlim(min(dist),max(dist))
    bathplot.set_ylim(min(deepthick),0)

    show()
