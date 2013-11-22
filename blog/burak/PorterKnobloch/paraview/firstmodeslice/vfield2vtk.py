#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np
import twomode

x = np.array([[0.75, 0, 0.1, 0.1],
			  [0.5, 0, 0.5, 0.5], 
			  [0.1, 0, 0.75, 0.75]], float)

T = twomode.generator()
			  
t = np.dot(T, x.transpose()).transpose()

npoints = np.size(x,0);

print "# vtk DataFile Version 3.0"
print "Vector field to define the tangents"
print "ASCII"

print "\nDATASET UNSTRUCTURED_GRID"
print "POINTS %i float" % npoints

for i in range(npoints):
	print x[i,0], x[i,1], x[i,2]

print "\nPOINT_DATA %i" % npoints
print "VECTORS tangentfield float "

for i in range(npoints):
	print t[i,0], t[i,1], t[i,2]
