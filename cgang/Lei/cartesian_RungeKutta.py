from matplotlib import pylab
from mpl_toolkits.mplot3d import axes3d
import numpy
import math

# time increment
dt = 0.00001

a1 = -0.4
a2 = -6
b1 = 1.6
b2 = -0.1
c1 = 1
c2 = 1
m1 = -0.337
m2 = 0.27
e = 1.219

def f_dx1(x1, y1, x2, y2):
    global a1, a2, b1, b2, c1, c2, m1, m2, e
    return x1*x2+y1*y2+x1*(m1+a1*(x1**2+y1**2)+b1*(x2**2+y2**2))

def f_dy1(x1, y1, x2, y2):
    global a1, a2, b1, b2, c1, c2, m1, m2, e
    return x1*y2-y1*x2+y1*(m1+a1*(x1**2+y1**2)+b1*(x2**2+y2**2))

def f_dx2(x1, y1, x2, y2):
    global a1, a2, b1, b2, c1, c2, m1, m2, e
    return c2*(x1**2-y1**2)+x2*(m2+a2*(x1**2+y1**2)+b2*(x2**2+y2**2))

def f_dy2(x1, y1, x2, y2):
    global a1, a2, b1, b2, c1, c2, m1, m2, e
    return e+c2*2*x1*y1+y2*(m2+a2*(x1**2+y1**2)+b2*(x2**2+y2**2))

	
#	foutrh order Runge-Kutta routine	
def rk4(x1,y1,x2,y2):
    global dt
    k1x1 = dt * f_dx1(x1,y1,x2,y2)
    k1y1 = dt * f_dy1(x1,y1,x2,y2)
    k1x2 = dt * f_dx2(x1,y1,x2,y2)
    k1y2 = dt * f_dy2(x1,y1,x2,y2)


    k2x1 = dt * f_dx1(x1 + k1x1/2.0, y1 + k1y1/2.0, x2 + k1x2/2.0, y2 + k1y2/2.0)
    k2y1 = dt * f_dy1(x1 + k1x1/2.0, y1 + k1y1/2.0, x2 + k1x2/2.0, y2 + k1y2/2.0)
    k2x2 = dt * f_dx2(x1 + k1x1/2.0, y1 + k1y1/2.0, x2 + k1x2/2.0, y2 + k1y2/2.0)
    k2y2 = dt * f_dy2(x1 + k1x1/2.0, y1 + k1y1/2.0, x2 + k1x2/2.0, y2 + k1y2/2.0)


    k3x1 = dt * f_dx1(x1 + k2x1/2.0, y1 + k2y1/2.0, x2 + k2x2/2.0, y2 + k2y2/2.0)
    k3y1 = dt * f_dy1(x1 + k2x1/2.0, y1 + k2y1/2.0, x2 + k2x2/2.0, y2 + k2y2/2.0)
    k3x2 = dt * f_dx2(x1 + k2x1/2.0, y1 + k2y1/2.0, x2 + k2x2/2.0, y2 + k2y2/2.0)
    k3y2 = dt * f_dy2(x1 + k2x1/2.0, y1 + k2y1/2.0, x2 + k2x2/2.0, y2 + k2y2/2.0)


    k4x1 = dt * f_dx1(x1 + k3x1, y1 + k3y1, x2 + k3x2, y2 + k3y2)
    k4y1 = dt * f_dy1(x1 + k3x1, y1 + k3y1, x2 + k3x2, y2 + k3y2)
    k4x2 = dt * f_dx2(x1 + k3x1, y1 + k3y1, x2 + k3x2, y2 + k3y2)
    k4y2 = dt * f_dy2(x1 + k3x1, y1 + k3y1, x2 + k3x2, y2 + k3y2)

    x1 = x1 + k1x1/6.0 + k2x1/3.0 + k3x1/3.0 +k4x1/6.0
    y1 = y1 + k1y1/6.0 + k2y1/3.0 + k3y1/3.0 +k4y1/6.0
    x2 = x2 + k1x2/6.0 + k2x2/3.0 + k3x2/3.0 +k4x2/6.0
    y2 = y2 + k1y2/6.0 + k2y2/3.0 + k3y2/3.0 +k4y2/6.0

    return [x1, y1, x2, y2]


all_x1=[]
all_y1=[]
all_x2=[]
all_y2=[]

# set initial values
i = 0
while i<1:
    x1=0.1+i
    y1=0.2+i
    x2=0.3+i
    y2=0.4+i
	
# list of points for one orbit
    list_x1 = [x1]
    list_y1 = [y1]
    list_x2 = [x2]
    list_y2 = [y2]

    t = 1
    while t < 100:
        [x1, y1, x2, y2] = rk4(list_x1[t-1], list_y1[t-1], list_x2[t-1], list_y2[t-1])
        print([x1,y1,x2,y2])
        list_x1.append(x1)
        list_y1.append(y1)
        list_x2.append(x2)
        list_y2.append(y2)
        t = t + 1
    all_x1=all_x1+list_x1
    all_y1=all_y1+list_y1
    all_x2=all_x2+list_x2
    all_y2=all_y2+list_y2
    i+=1
	





# draw graph using
fig = pylab.figure()
ax = axes3d.Axes3D(fig)
ax.scatter3D(all_x1, all_x2, all_y2, color='black', s = 0.1)
#ax.scatter3D(list_x1,list_y1, list_x2, s = 0.1)

ax.set_xlabel('X1')
ax.set_ylabel('Y1')
ax.set_zlabel('X2')
pylab.show()
