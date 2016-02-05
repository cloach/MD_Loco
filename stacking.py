#!/usr/bin/python

import sys
import math

argument_list = sys.argv
element = argument_list[1]
stacking = argument_list[2]
lattice_constant = float(argument_list[3])

interplanar_spacing = lattice_constant/math.sqrt(3)
layers = len(stacking)
height = layers*interplanar_spacing
width = lattice_constant/math.sqrt(2)
 
print '%Block Lattice_ABC'
print width, width, height
print '90 90 60'
print '%EndBlock Lattice_ABC'
print ''
print '%Block Positions_Frac'
for i in range (0, len(stacking)):
    if stacking[i]=='a':
	    print element, float(0), float(0), float(i)/float(layers)
    elif stacking[i]=='b':
	    print element, float(1)/float(3), float(1)/float(3), float(i)/float(layers)
    elif stacking[i]=='c':
	    print element, float(2)/float(3), float(2)/float(3), float(i)/float(layers)
print 'EndBlock Positions_Frac'

