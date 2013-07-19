import math
a1 = 1.
a2 = 1.5
b1 = 3.
b2 = 2.5
c1 = 3.
c2 = 3.5
m1 = -1.
m2 = -4.
e = 0.1

r1 = 0.467095
r2 = 0.2146
psi = e/(c2*r1**2/r2+2*c1*r2) 

a=[[0,0,0],[0,0,0],[0,0,0]]
a[0][0] = m1+2*a1*r1**2+b1*r2**2+c1*r2*math.cos(psi)
a[0][1] = 2*b1*r1*r2+c1*r1*math.cos(psi)
a[0][2] = -c1*r1*r2*math.sin(psi)
a[1][0] = 2*c2*r1*math.cos(psi)
a[1][1] = m2+a1*r1**2+3*b2*r2**2
a[1][2] = -c2*r1**2*math.sin(psi)
a[2][0] = -c2*2*r1/r2*math.sin(psi)
a[2][1] = c2*r1**2/(r2**2)*math.sin(psi)-2*c1*math.sin(psi)
a[2][2] = -(c2*r1**2/r2+2*c1*r2)*math.cos(psi)
for i in range(3):
    for j in range(3):
	    print a[i][j]
