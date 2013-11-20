#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np
import twomode

x = np.loadtxt('data/solutiononslice.dat')

npoints = np.size(x,0);
npoints = 50000

print "# vtk DataFile Version 3.0"
print "Line data to hold the curve in 3D"
print "ASCII"

print "\nDATASET POLYDATA"
print "POINTS %i float" % npoints

for i in range(npoints):
	print x[i,1], x[i,3], x[i,4]

print "LINES 1 %i" % (npoints+1)

for i in range(npoints+1):
	j = (i - 1) % (npoints + 1)
	print(j),
