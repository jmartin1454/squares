# squares

To install:

1.  git clone
2.  cd squares
3.  git submodule init
4.  git submodule update

Then run squarespeed.py.

What the code does:
- sets up an array of coils inscribed on a cube
- sets up an array of fluxgates distributed throughout a smaller cube
- calculates the matrix in B=M*I by setting each I in turn and measuring B
- inverts the matrix, carefully excising the zero mode.
- selects a target field B_target (uses the Pis library)
- calculates I_set = M^{-1} B_target and sets those currents on the cube
- makes graphs of the resultant true field and graphs them in a few views (presently along axes), comparing with the target field.
