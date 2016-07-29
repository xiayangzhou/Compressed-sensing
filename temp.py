 -*- coding: iso-8859-1 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from cvxopt import matrix
from AGKM2 import AGKM
from cvxopt import normal
import numpy as np
A = matrix([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], (2,3))
B = matrix ([5,6,7],(1,3))
print(A)

#A.size

a, b, c = 500, 100, 10
#F ,W= np.normal(b,c),np.normal(c,a)

W1=np.matrix([[1,0,0.,0],[0,1,0,0]])
F1=np.matrix([[1,0],[.5,.5],[0.25,0.75],[0,1]])
Y = np.dot(F1,W1)
X = matrix(Y)

Q=AGKM(X,0.5,0.0001)

