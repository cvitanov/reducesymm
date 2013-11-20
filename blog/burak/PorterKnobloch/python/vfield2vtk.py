#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np

x = np.array([[0,0,0], 
			  [0.4399655797367152, -0.38626706847930564, 0.0702043939917171]], float)
			  
t = np.array([[0, 0.70710678, -0.70710678], 
			  [0.02884567, 0.99957642, -0.00386173]], float)

npoints = np.size(x,0);

print "# vtk DataFile Version 3.0"
print "Vector field to define the tangents"
print "ASCII"

print "\nDATASET UNSTRUCTURED_GRID"
print "POINTS %i float" % npoints

for i in range(npoints):
	print x[i,0], x[i,1], x[i,2]

print "\nPOINT_DATA %i" % npoints
print "VECTORS unstabledirs float "

for i in range(npoints):
	print t[i,0], t[i,1], t[i,2]
