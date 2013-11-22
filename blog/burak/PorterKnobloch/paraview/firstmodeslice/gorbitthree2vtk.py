#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np
import twomode

x = np.array([0.1, 0, 0.75, 0.75], float)
theta = np.arange(0, 2*np.pi, 0.001)

gx	= np.zeros((np.size(theta,0),4))			  

for i in range(0, np.size(theta,0), 1): 
	gx[i,:] = np.dot(twomode.LieElement(theta[i]), x)

npoints = np.size(gx,0);

print "# vtk DataFile Version 3.0"
print "Line data to hold the curve in 3D"
print "ASCII"

print "\nDATASET POLYDATA"
print "POINTS %i float" % npoints

for i in range(npoints):
	print gx[i,0], gx[i,1], gx[i,2]

print "LINES 1 %i" % (npoints+1)

for i in range(npoints+1):
	j = (i - 1) % (npoints + 1)
	print(j),
