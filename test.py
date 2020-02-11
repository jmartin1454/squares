#!/usr/bin/python

# Does many of the manipulations described in:
# https://www.cs.bgu.ac.il/~na181/wiki.files/SVD_application_paper[1].pdf
# but
# using the python libraries and methods described in
# https://machinelearningmastery.com/singular-value-decomposition-for-machine-learning/

import numpy as np

from numpy import array
from numpy.linalg import svd
from numpy.linalg import pinv
from numpy import zeros
from numpy import diag
# define matrix
A = array([
	[1., 2.],
	[2., 4.],
	[3., 6.],
	[4., 8.0001]])
print("A",A)
lilb=array([2,
            4,
            6.000,
            8])
# calculate svd
U, s, VT = svd(A)
print("U",U)
print("s",s)
print("VT",VT)
# reciprocals of s
d = 1.0 / s
# create m x n D matrix
D = zeros(A.shape)
# populate D with n x n diagonal matrix
D[:A.shape[1], :A.shape[1]] = diag(d)
print("D",D)
# calculate pseudoinverse
B = VT.T.dot(D.T).dot(U.T)
print("B",B)
print("Bb",B.dot(lilb))
print(pinv(A))
print(pinv(A).dot(lilb))
# fix by cutting on cond #
print("fixed",pinv(A,.0001))
print(pinv(A,.0001).dot(lilb))
print("dodo")
# do same for svd
# select
n_elements = 1 # fix by cutting to this number of eigenvalues
D = D[:, :n_elements]
VT = VT[:n_elements, :]
print("newD",D)
print("newVT",VT)
B = VT.T.dot(D.T).dot(U.T)
print("newB",B)
print("newBb",B.dot(lilb))
