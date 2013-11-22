#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np
import twomode

center = np.array([0.5, 0, 0.25], float)
basis1 = np.array([0.5,  0, 0], float)/1
basis2 = np.array([0,0,0.8], float)/1

x = np.array([center-basis1-basis2, center+basis1-basis2, center+basis1+basis2, center-basis1+basis2], float)

npoints = np.size(x,0);

print "# vtk DataFile Version 3.0"
print "Line data to hold the curve in 3D"
print "ASCII"

print "\nDATASET POLYDATA"
print "POINTS %i float" % npoints

for i in range(npoints):
	print x[i,0], x[i,1], x[i,2]

print "POLYGONS 1 %i" % (npoints+1)

for i in range(npoints+1):
	j = (i - 1) % (npoints + 1)
	print(j),
