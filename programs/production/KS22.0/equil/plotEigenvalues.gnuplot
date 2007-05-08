set term postscript eps enhanced color  		 # enhanced PostScript, essentially PostScript
 		 			 # with bounding boxes
set out ' L22-eqvaEigenvalues.eps'

set xrange [-2.2:0.3]    
set grid

set xlabel '{/Symbol m}_i' 
set ylabel '{/Symbol n}_i' 
#set lmargin 10
#set label 1 '{/Symbol n}_i' at graph -0.2, graph 0.5
set key top left 

set pointsize 2

plot "1w/Jdiag.dat" using 1:2 title 'E_1',"2w/Jdiag.dat" using 1:2 title 'E_2',"3w/Jdiag.dat" using 1:2 title 'E_3' 
