#
#symbdyn2tex.py
#
"""
	Python script to read data/symbdyn.dat and convert it to a latex table
	Used approach: http://stackoverflow.com/questions/7111690/python-read-formatted-string
"""

import numpy as np

xpo = np.loadtxt('data/xrpo.dat')
itineraries = np.loadtxt('data/itineraries.dat', dtype="str")
periods = np.loadtxt('data/periods.dat')
position = np.loadtxt('data/position.dat')
group = np.loadtxt('data/group.dat')
tofpo = np.loadtxt('data/tofrpo.dat')

lines = []

for i in range(1,int(np.max(group))+1):
	
	#indices corresponding to the members of ith periodic orbit:
	gindices = np.argwhere(group == i)
	gindices = gindices.reshape(np.size(gindices))
	print "gindices:"
	print gindices
	#Position of first element of the cycle:
	p1 = np.argwhere(position[gindices] == 1)
	p1 = p1.reshape(np.size(p1))
	print "p1:"
	print p1
	#Position of the first element:
	i1 = gindices[p1]
	i1 = i1.reshape(np.size(i1))
	print "i1:"
	print i1
	x0 = xpo[i1,:]
	x0 = x0.reshape(np.size(x0))
	print "x0:"
	print x0
	itinerary = itineraries[i1]
	print itinerary[0]
	
	#Period of the ith cycle:
	T = periods[i-1]
	print "Period ="
	print T
	
	lines.append([itinerary[0], x0[0], x0[1], x0[2], x0[3], T])

print lines

#raw_input("Press Enter to continue...")

f = open("data/twomoderpos.tex", "w")

f.write("\\begin{table}\n")
f.write("\t\\begin{tabular}{c|c|c}\n")

f.write("\tItinerary & $(x_{1,RPO}, y_{1,RPO}, x_{2,RPO}, y_{2,RPO})$ & Period \\\\ \n")

for i in range(len(lines)):
	
	f.write("\t\\hline\n")
	f.write("\t%s & " % lines[i][0])
	f.write("(%5.7f, " % float(lines[i][1])) 
	f.write("%5.1f, " % abs(float(lines[i][2]))) 
	f.write("%5.7f, " % float(lines[i][3])) 
	f.write("%5.7f) & " % float(lines[i][4])) 
	f.write("%5.7f \\\\ \n " % float(lines[i][5])) 
	
f.write("\t\\end{tabular}\n")
f.write("\t\\caption{\\rpo s of the \\twoMode\\ system.}\n")
f.write("\t\\label{tab:twomoderpos}\n")
f.write("\\end{table}")

f.close()
