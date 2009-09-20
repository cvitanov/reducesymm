# $Author: predrag $ $Date: 2006-07-06 19:05:16 -0400 (Thu, 06 Jul 2006) $
# siminos/figSrc/gnu/lyapSpec.gnu

# PC 2009-09-13: this plots Ruslan's KS L=22 Lyapunov exponents
#                lyapSpec.eps, used in siminos/blog/blog.tex

reset
#set term x11
set term postscript eps enhanced color solid "Times-Roman" 22
# set terminal fig
# set output "lyapSpec.fig"

set title "KS L=22 Lyapunov spectrum" font "Times-Roman,24"
set output "lyapSpecRscld.eps" # in j 2 \pi/L scale

# color size 20 10 metric solid fontsize 16 textnormal thickness 2k depth 30
# {pointsmax <max_points>}

# set label "f_0" at 0.5,0.3
set size 1,1
set border 15
set xrange [0:2.1]
set yrange [-12.5:0.95]
set xlabel "2{/Symbol p}j/L"
# PC: this is postscriptese for \lambda_j:
set ylabel "{/Symbol l}_{j}"
set pointsize 2
# PC: points -> circles is obtained by point_type 6
#            -> + is obtained by point_type 1

g1(x)=x**2-x**4
#PC example:  g1(x)=exp(-(x-a)*(x-a)/g/g/2.0)/sqrt(pi)/g

L=22  # PC: plot in `extensive,' discretization independent scale
tildeL = L/(2*pi) 
#PC: in lyapSpecRscld.dat j runs over pairs of eigenvalues

plot g1(x) notitle, \
"lyapSpecRscld.dat" using ($1/tildeL):($2) notitle with points pointtype 1, \
"lyapSpecEQ3.dat" using ($1/tildeL):($2) notitle with points pointtype 1


# set style line 1 lt 1 lw 3
# set style line 2 lt 1 lw 1
# set arrow from a,0.35 to a,0.0 nohead linestyle 1
# set arrow from c,p1   to 14,p1 linestyle 2
# set label "contribution from A" at 14,p1

set xtics 
set ytics
# set xtics nomirror ("E_A" a, "E_x" c, "E_B" b)
# set ytics nomirror ("" 0)
