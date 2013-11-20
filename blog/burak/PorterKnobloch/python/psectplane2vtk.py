#
# plot2vtk.py
#
"""
Generates a vtk file to hold 3D vector field as a vtk file in
UNSTRUCTURED_GRID format
"""

import numpy as np
import twomode

eq1 = np.array([0.4399655797367152, -0.38626706847930564, 0.0702043939917171], float)
vu = np.array([0.02884567,  0.99957642, -0.00386173], float)/1
vaux = np.array([0,0,1], float)/2

x = np.array([eq1-vu-vaux, eq1+vu-vaux, eq1+vu+vaux, eq1-vu+vaux], float)

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
