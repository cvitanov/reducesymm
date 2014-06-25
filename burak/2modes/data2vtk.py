"""
Read the data created by main and produce vtk files to be imported from Paraview
"""
import numpy as np
import sys

txphisol = np.loadtxt('data/txphisol.dat')
npoints = np.size(txphisol,0);

orig_stdout = sys.stdout

f = file('vtk/flow.vtk', 'w')
sys.stdout = f

print "# vtk DataFile Version 3.0"
print "Line data to hold the curve in 3D"
print "ASCII"

print "\nDATASET POLYDATA"
print "POINTS %i float" % npoints

for i in range(npoints):
	print txphisol[i,1], txphisol[i,3], txphisol[i,4]

print "LINES 1 %i" % (npoints+1)

for i in range(npoints+1):
	j = (i - 1) % (npoints + 1)
	print(j),

sys.stdout = orig_stdout
f.close()

ps = np.loadtxt('data/ps.dat')
npoints = np.size(ps,0)


f = file('vtk/ps.vtk', 'w')
sys.stdout = f

print "# vtk DataFile Version 3.0"
print "Line data to hold the curve in 3D"
print "ASCII"

print "\nDATASET UNSTRUCTURED_GRID"
print "POINTS %i float" % npoints

for i in range(npoints):
	print ps[i,1], ps[i,3], ps[i,4]

sys.stdout = orig_stdout
f.close()

reqv = np.loadtxt('data/reqv3d.dat')
unstabledir = np.loadtxt('data/unstabledir3d.dat')
unstabledir2 = np.loadtxt('data/unstabledir23d.dat')

f = file('vtk/unstabledir.vtk', 'w')
sys.stdout = f

print "# vtk DataFile Version 3.0"
print "Vector field to define the tangents"
print "ASCII"

print "\nDATASET UNSTRUCTURED_GRID"
print "POINTS %i float" % 2

print reqv[0], reqv[1], reqv[2]
print reqv[0], reqv[1], reqv[2]

print "\nPOINT_DATA %i" % 2
print "VECTORS unstabledirs float "

print unstabledir[0], unstabledir[1], unstabledir[2]
print unstabledir2[0], unstabledir2[1], unstabledir2[2]

sys.stdout = orig_stdout
f.close()

center = reqv
basis1 = np.array([0,  0, 1], float)/4
basis2 = unstabledir*2.5

#x = np.array([center-basis1-basis2, center+basis1-basis2, center+basis1+basis2, center-basis1+basis2], float)
x = np.array([center-basis1, center+basis1, center+basis1+basis2, center-basis1+basis2], float)
npoints = np.size(x,0);

f = file('vtk/psectplane.vtk', 'w')
sys.stdout = f

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

sys.stdout = orig_stdout
f.close()
