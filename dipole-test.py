#!/usr/bin/python3

from dipole import dipole

# magnetic dipole moment in A*m**2
mx=0
my=0
mz=1

# location of the dipole
x0=0
y0=0
z0=0

# initialize the dipole with this magnetic moment and location
d=dipole(x0,y0,z0,mx,my,mz)

# what's the field produced by this dipole at some location?
x=1
y=0
z=0

print(d.bx(x,y,z),d.by(x,y,z),d.bz(x,y,z),d.b(x,y,z))
