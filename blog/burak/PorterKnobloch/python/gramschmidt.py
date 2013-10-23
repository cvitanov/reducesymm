#
# movingframes3dplotter.py
#
"""
Convert solutions to Gram-Schmidt coordinates
"""

from numpy import loadtxt
import numpy as np
import twomode

t, x1hat, x2hat, y1hat, y2hat, phi = loadtxt('data/solutiononslice.dat', unpack=True)

x = np.array([x1hat, x2hat, y1hat, y2hat], float)
xgs = twomode.ssp2gramschmidt(x)
xgs = xgs.transpose()

#Print it:
for t1, x1 in zip(t,xgs):
	print t1, x1[0], x1[1], x1[2]
