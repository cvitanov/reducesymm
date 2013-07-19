(* ::Package:: *)

(*Given a set of parameters, find equilibrium points and corresponding eigenvalues*)

a1 = -0.245
a2 = -3.44
b1 = 1.326
b2 = -0.47
c1 = 1
c2 = -1
m1 = -0.14
m2 = 1.175
e = 0.855

sol={x,y} /. NSolve[{(m1*x+a1*x^3+b1*x*y^2)^2/(c1^2*x^2*y^2)+(y^2*e^2)/(c2*x^2+2*c1*y^2)^2==1,(m2*y+a2*x^2*y+b2*y^3)^2/(c2^2*x^4)+(y^2*e^2)/(c2*x^2+2*c1*y^2)^2==1, x>=0, y>=0},{x,y},Reals]



(*first equilibrium*)
r1 = sol[[1,1]]
r2 = sol[[1,2]]
psi = e/(c2*r1^2/r2+2*c1*r2) 

NSolve[{r1/x*r2*y+r1*x*(m1+a1*r1^2+b1*r2^2)==0,c2*r1*x^2+r2*y*(m2+a2*r1^2+b2*r2^2)==0},{x,y}]

a={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
a[[1,1]] = m1+2*a1*r1^2+b1*r2^2+c1*r2*Cos[psi]
a[[1,2]] = 2*b1*r1*r2+c1*r1*Cos[psi]
a[[1,3]] = -c1*r1*r2*Sin[psi]
a[[2,1]] = 2*c2*r1*Cos[psi]
a[[2,2]] = m2+a1*r1^2+3*b2*r2^2
a[[2,3]] = -c2*r1^2*Sin[psi]
a[[3,1]] = -c2*2*r1/r2*Sin[psi]
a[[3,2]] = c2*r1^2/(r2^2)*Sin[psi]-2*c1*Sin[psi]
a[[3,3]] = -(c2*r1^2/r2+2*c1*r2)*Cos[psi]

eigs[1] = Eigenvalues[a]


(*second equilibrium*)
r1 = solp[[2,1]]
r2 = solp[[2,2]]
psi = e/(c2*r1^2/r2+2*c1*r2) 

a={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
a[[1,1]] = m1+2*a1*r1^2+b1*r2^2+c1*r2*Cos[psi]
a[[1,2]] = 2*b1*r1*r2+c1*r1*Cos[psi]
a[[1,3]] = -c1*r1*r2*Sin[psi]
a[[2,1]] = 2*c2*r1*Cos[psi]
a[[2,2]] = m2+a1*r1^2+3*b2*r2^2
a[[2,3]] = -c2*r1^2*Sin[psi]
a[[3,1]] = -c2*2*r1/r2*Sin[psi]
a[[3,2]] = c2*r1^2/(r2^2)*Sin[psi]-2*c1*Sin[psi]
a[[3,3]] = -(c2*r1^2/r2+2*c1*r2)*Cos[psi]

eigs[2] = Eigenvalues[a]