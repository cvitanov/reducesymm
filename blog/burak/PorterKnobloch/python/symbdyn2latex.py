#
#symbdyn2tex.py
#
"""
	Python script to read data/symbdyn.dat and convert it to a latex table
	Used approach: http://stackoverflow.com/questions/7111690/python-read-formatted-string
"""

f = open("data/symbdyn.dat", "r")

lines = []

for line in f:
	lines.append(line.split("\t"))
	
f.close()

#print lines

f = open("data/twomodesymbdynBB.tex", "w")

f.write("\\begin{table}\n")
f.write("\t\\begin{tabular}{c|c}\n")

f.write("\tItinerary & $(x_{1,RPO}, y_{1,RPO}, x_{2,RPO}, y_{2,RPO})$ \\\\ \n")

for i in range(len(lines)):
	
	f.write("\t\\hline\n")
	f.write("\t%s & " % lines[i][0])
	f.write("(%5.7f, " % float(lines[i][1])) 
	f.write("%5.7f, " % float(lines[i][2])) 
	f.write("%5.7f, " % float(lines[i][3])) 
	f.write("%5.7f) \\\\ \n " % float(lines[i][4])) 
	
f.write("\t\\end{tabular}\n")
f.write("\t\\caption{Symbolic dynamics of \\twoMode\\ system.}\n")
f.write("\t\\label{tab:symbdyn}\n")
f.write("\\end{table}")

f.close()
