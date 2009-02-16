#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import cos, sin, atan2, radians, degrees, sqrt, log, pi

x = 0
def func1():
    global x
    print x    
    x = 22

def func2():
   global x
   print x

func1()
func2()
