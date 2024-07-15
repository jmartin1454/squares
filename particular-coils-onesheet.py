#!/usr/bin/env python3

# Fri May 24 10:11:45 CDT 2019 Jeff added this line.

# Tue Feb 11 13:43:43 CST 2020 Jeff taking original patch.py and
# updating to solve the zero mode issue.  Will now update to use the
# patchlib submodule.

# Fri Feb 14 11:45:17 CST 2020 Jeff speeding up code by improving the
# sensorarray class to use numpy structures.

# Sat Aug 1 12:40:53 CDT 2020 Jeff added command line options and
# improved graphs


from scipy.constants import mu_0, pi
import numpy as np
from patchlib.patch import *
from Pis.Pislib import *
from dipole import *

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-s", "--nsensors", dest="nsensors", default=3,
                  help="ns where total sensor axes is s = 3*ns^3")

parser.add_option("-l", "--ell", dest="l", default=2,
                  help="l for spherical harmonic")

parser.add_option("-m", "--em", dest="m", default=0,
                  help="m for spherical harmonic")

parser.add_option("-M", "--matrices", dest="matrices", default=False,
                  action="store_true",
                  help="show matrices")

parser.add_option("-d", "--dipole", dest="dipole", default=False,
                  action="store_true",
                  help="use dipole field")

parser.add_option("-t", "--traces", dest="traces", default=False,
                  action="store_true",
                  help="show 3D view of coils and sensors")

parser.add_option("-r", "--residuals", dest="residuals", default=False,
                  action="store_true",
                  help="show residuals")

parser.add_option("-z", "--zoom", dest="zoom", default=False,
                  action="store_true",
                  help="zoom to ROI")

parser.add_option("-a", "--axes", dest="axes", default=False,
                  action="store_true",
                  help="make graphs along axes")

parser.add_option("-i", "--incells", dest="incells", default=False,
                  action="store_true",
                  help="ROI for statistics is in EDM cells")

#d=dipole(1.2,0,0,0,0,100000)  # dipole1
d=dipole(0,0,1.2,0,0,1)  # dipole2
#d=dipole(0,0,1.2,1,0,0)  # dipole3

(options,args)=parser.parse_args()

l=int(options.l)
m=int(options.m)
sp=scalarpotential(l,m)
print("Sigma in spherical coordinates is %s"%sp.Sigma_spherical)
print("Sigma in cartesian coordinates is %s"%sp.Sigma)

print("Pix is %s"%sp.Pix)
print("Piy is %s"%sp.Piy)
print("Piz is %s"%sp.Piz)

if(options.dipole):
    bxtarget=d.bx
    bytarget=d.by
    bztarget=d.bz
else:
    bxtarget=sp.fPix
    bytarget=sp.fPiy
    bztarget=sp.fPiz

# Setup our coilset

myset=coilset()

# positions of faces (positive values -- faces will be at plus and
# minus of these)
a=2.2 # approximate position of layer 6 of MSR
xface=a/2 # m
yface=a/2 # m
zface=a/2 # m

# set up the rear wall
y1=250*.001
z1=260*.001
x1=xface
point1=(x1,y1,z1)
y2=y1
z2=600*0.001  # Jeff changed to fit on 4' sheet
x2=xface
point2=(x2,y2,z2)
y3=600*.001 # guess
z3=z2
x3=xface
point3=(x3,y3,z3)
y4=y3
z4=z1
x4=xface
point4=(x4,y4,z4)
points_ur=(point1,point2,point3,point4)
points_ur=np.array(points_ur)
myset.add_coil(points_ur)

# Now add mirror images of these
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
points_ul=(point1,point2,point3,point4)
points_ul=np.array(points_ul)
myset.add_coil(points_ul)

point1=(x1,-y1,-z1)
point2=(x2,-y2,-z2)
point3=(x3,-y3,-z3)
point4=(x4,-y4,-z4)
points_ll=(point1,point2,point3,point4)
points_ll=np.array(points_ll)
myset.add_coil(points_ll)

point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
points_lr=(point1,point2,point3,point4)
points_lr=np.array(points_lr)
myset.add_coil(points_lr)

# now the central coil
x1=xface
#y1=530/2*.001
y1=190*.001
z1=400/2*.001
point1=(x1,y1,z1)
point2=(x1,y1,-z1)
point3=(x1,-y1,-z1)
point4=(x1,-y1,z1)
points_c=(point1,point2,point3,point4)
points_c=np.array(points_c)
myset.add_coil(points_c)

# now the right side coil
x1=xface
y1=250*0.001
z1=400/2*.001
point1=(x1,y1,z1)
x2=xface
#y2=y1+420*.001 # guess
y2=600*.001 # guess
z2=z1
point2=(x2,y2,z2)
x3=xface
y3=y2
z3=-z2
point3=(x3,y3,z3)
x4=xface
y4=y1
z4=-z1
point4=(x4,y4,z4)
points_mr=(point1,point2,point3,point4)
print('points_mr',points_mr)
points_mr=np.array(points_mr)
myset.add_coil(points_mr)

# now the left side coil -- reflect and wind in same direction
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
points_ml=(point1,point2,point3,point4)
points_ml=np.array(points_ml)
myset.add_coil(points_ml)

# now the upper central coil
x1=xface
y1=190*0.001
z1=260*.001
point1=(x1,y1,z1)
x2=xface
y2=-190*.001 # guess
z2=z1
point2=(x2,y2,z2)
x3=xface
y3=y2
z3=600*.001
point3=(x3,y3,z3)
x4=xface
y4=y1
z4=z3
point4=(x4,y4,z4)
points_uc=(point1,point2,point3,point4)
print('points_uc',points_uc)
points_uc=np.array(points_uc)
myset.add_coil(points_uc)

# now the lower central coil -- reflect and wind in same direction
point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
points_lc=(point1,point2,point3,point4)
points_lc=np.array(points_lc)
myset.add_coil(points_lc)

# now reflect them all to the other face: xface -> -xface
def reflect_x(points):
    newpoints=np.copy(points)
    newpoints[:,0]=-newpoints[:,0]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints
    
oside_ur=reflect_x(points_ur)
myset.add_coil(oside_ur)
oside_ul=reflect_x(points_ul)
myset.add_coil(oside_ul)
oside_ll=reflect_x(points_ll)
myset.add_coil(oside_ll)
oside_lr=reflect_x(points_lr)
myset.add_coil(oside_lr)
oside_c=reflect_x(points_c)
myset.add_coil(oside_c)
oside_ml=reflect_x(points_ml)
myset.add_coil(oside_ml)
oside_mr=reflect_x(points_mr)
myset.add_coil(oside_mr)
oside_uc=reflect_x(points_uc)
myset.add_coil(oside_uc)
oside_lc=reflect_x(points_lc)
myset.add_coil(oside_lc)

# Phew -- now onto the sides

z1=(400/2+60)*.001
x1=(105/2+60+255+60)*.001
y1=-yface
point1=(x1,y1,z1)
z2=z1+420*.001 # guess
x2=x1
y2=-yface
point2=(x2,y2,z2)
z3=z2
x3=x2+500*.001 # guess
y3=-yface
point3=(x3,y3,z3)
z4=z1
x4=x3
y4=-yface
point4=(x4,y4,z4)
side_ur=(point1,point2,point3,point4)
side_ur=np.array(side_ur)
myset.add_coil(side_ur)

# now reflect around
point1=(-x1,y1,z1)
point2=(-x4,y4,z4)
point3=(-x3,y3,z3)
point4=(-x2,y2,z2)
side_ul=np.array((point1,point2,point3,point4))
myset.add_coil(side_ul)

point1=(-x1,y1,-z1)
point2=(-x2,y2,-z2)
point3=(-x3,y3,-z3)
point4=(-x4,y4,-z4)
side_ll=np.array((point1,point2,point3,point4))
myset.add_coil(side_ll)

point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
side_lr=np.array((point1,point2,point3,point4))
myset.add_coil(side_lr)

# central coil
z1=400/2*.001
y1=-yface
x1=(105/2+60+255)*.001
point1=(x1,y1,z1)
point2=(x1,y1,-z1)
point3=(-x1,y1,-z1)
point4=(-x1,y1,z1)
side_c=np.array((point1,point2,point3,point4))
myset.add_coil(side_c)

# middle right coil
x1=(105/2+60+255+60)*.001
y1=-yface
z1=400/2*.001
point1=(x1,y1,z1)
x2=x1+500*.001 # same guess as above
y2=-yface
z2=z1
point2=(x2,y2,z2)
point3=(x2,y2,-z2)
point4=(x1,y1,-z1)
side_mr=np.array((point1,point2,point3,point4))
myset.add_coil(side_mr)

# reflect it to middle left coil
point1=(-x1,y1,z1)
point2=(-x1,y1,-z1)
point3=(-x2,y2,-z2)
point4=(-x2,y2,z2)
side_ml=np.array((point1,point2,point3,point4))
myset.add_coil(side_ml)

# middle top
z1=(400/2+60)*.001
x1=(105/2+60+255)*.001
y1=-yface
point1=(x1,y1,z1)
z2=z1
x2=-x1
y2=-yface
point2=(x2,y2,z2)
z3=z2+420*.001 # same guess as above
x3=x2
y3=-yface
point3=(x3,y3,z3)
z4=z3
x4=x1
y4=-yface
point4=(x4,y4,z4)
side_mt=np.array((point1,point2,point3,point4))
myset.add_coil(side_mt)

# mirror to middle bottom
point1=(x1,y1,-z1)
point2=(x4,y4,-z4)
point3=(x3,y3,-z3)
point4=(x2,y2,-z2)
side_mb=np.array((point1,point2,point3,point4))
myset.add_coil(side_mb)

# now reflect them all to the other face: -yface -> yface
def reflect_y(points):
    newpoints=np.copy(points)
    newpoints[:,1]=-newpoints[:,1]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints

oside_side_ur=reflect_y(side_ur)
oside_side_ul=reflect_y(side_ul)
oside_side_ll=reflect_y(side_ll)
oside_side_lr=reflect_y(side_lr)
oside_side_c=reflect_y(side_c)
oside_side_ml=reflect_y(side_ml)
oside_side_mr=reflect_y(side_mr)
oside_side_mt=reflect_y(side_mt)
oside_side_mb=reflect_y(side_mb)

myset.add_coil(oside_side_ur)
myset.add_coil(oside_side_ul)
myset.add_coil(oside_side_lr)
myset.add_coil(oside_side_ll)
myset.add_coil(oside_side_c)
myset.add_coil(oside_side_ml)
myset.add_coil(oside_side_mr)
myset.add_coil(oside_side_mt)
myset.add_coil(oside_side_mb)


# Double phew, now on to the top side

x1=(400/2+60)*.001 # picture frame of 400x400's separated by 60's
y1=(400/2+60)*.001
z1=zface
point1=(x1,y1,z1)
x2=x1
y2=y1+400*.001
z2=zface
point2=(x2,y2,z2)
x3=x2+400*.001
y3=y2
z3=zface
point3=(x3,y3,z3)
x4=x3
y4=y1
z4=zface
point4=(x4,y4,z4)
top_ur=(point1,point2,point3,point4)
top_ur=np.array(top_ur)
myset.add_coil(top_ur)

# now reflect around
point1=(-x1,y1,z1)
point2=(-x4,y4,z4)
point3=(-x3,y3,z3)
point4=(-x2,y2,z2)
top_ul=np.array((point1,point2,point3,point4))
myset.add_coil(top_ul)

point1=(-x1,-y1,z1)
point2=(-x2,-y2,z2)
point3=(-x3,-y3,z3)
point4=(-x4,-y4,z4)
top_ll=np.array((point1,point2,point3,point4))
myset.add_coil(top_ll)

point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
top_lr=np.array((point1,point2,point3,point4))
myset.add_coil(top_lr)

# central coil
z1=zface
y1=400/2*.001
x1=400/2*.001
point1=(x1,y1,z1)
point2=(x1,-y1,z1)
point3=(-x1,-y1,z1)
point4=(-x1,y1,z1)
top_c=np.array((point1,point2,point3,point4))
myset.add_coil(top_c)

# middle right coil
x1=(400/2+60)*.001
y1=400/2*.001
z1=zface
point1=(x1,y1,z1)
x2=x1+400*.001
y2=y1
z2=zface
point2=(x2,y2,z2)
point3=(x2,-y2,z2)
point4=(x1,-y1,z1)
top_mr=np.array((point1,point2,point3,point4))
myset.add_coil(top_mr)

# reflect it to middle left coil
point1=(-x1,y1,z1)
point2=(-x1,-y1,z1)
point3=(-x2,-y2,z2)
point4=(-x2,y2,z2)
top_ml=np.array((point1,point2,point3,point4))
myset.add_coil(top_ml)

# middle top
x1=400/2*.001
y1=(400/2+60)*.001
z1=zface
point1=(x1,y1,z1)
x2=-x1
y2=y1
z2=zface
point2=(x2,y2,z2)
x3=x2
y3=y2+400*.001
z3=zface
point3=(x3,y3,z3)
x4=x1
y4=y3
z4=zface
point4=(x4,y4,z4)
top_mt=np.array((point1,point2,point3,point4))
myset.add_coil(top_mt)

# mirror to middle bottom
point1=(x1,-y1,z1)
point2=(x4,-y4,z4)
point3=(x3,-y3,z3)
point4=(x2,-y2,z2)
top_mb=np.array((point1,point2,point3,point4))
myset.add_coil(top_mb)

# now reflect them all to the other face: zface -> -zface
def reflect_z(points):
    newpoints=np.copy(points)
    newpoints[:,2]=-newpoints[:,2]
    newpoints=np.flip(newpoints,0) # wind them in the opposite direction
    return newpoints

bott_ur=reflect_z(top_ur)
bott_ul=reflect_z(top_ul)
bott_ll=reflect_z(top_ll)
bott_lr=reflect_z(top_lr)
bott_c=reflect_z(top_c)
bott_ml=reflect_z(top_ml)
bott_mr=reflect_z(top_mr)
bott_mt=reflect_z(top_mt)
bott_mb=reflect_z(top_mb)

myset.add_coil(bott_ur)
myset.add_coil(bott_ul)
myset.add_coil(bott_lr)
myset.add_coil(bott_ll)
myset.add_coil(bott_c)
myset.add_coil(bott_ml)
myset.add_coil(bott_mr)
myset.add_coil(bott_mt)
myset.add_coil(bott_mb)


class sensor:
    def __init__(self,pos):
        self.pos = pos

class sensorarray:
    def __init__(self,xdim,ydim,zdim,corners):
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        #self.sensorgrid=np.mgrid[-a:a:xdim*1j,-a:a:ydim*1j,-a:a:zdim*1j]
        #print(self.sensorgrid)
        self.sensors = []
        if(xdim==1 and ydim==1 and zdim==1):
            pos = corners[0]+x/2+y/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+y/2+z
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2+x
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2+y
            self.sensors.append(sensor(pos))
        else:
            for i in range(xdim):
                for j in range(ydim):
                    for k in range(zdim):
                        pos = corners[0]+x*i/(xdim-1)+y*j/(ydim-1)+z*k/(zdim-1)
                        self.sensors.append(sensor(pos))
        self.numsensors = len(self.sensors)
    def draw_sensor(self,number,ax):
        x = self.sensors[number].pos[0]
        y = self.sensors[number].pos[1]
        z = self.sensors[number].pos[2]
        c = 'r'
        m = 'o'
        ax.scatter(x,y,z,c=c,marker=m)
    def draw_sensors(self,ax):
        for number in range(self.numsensors):
            self.draw_sensor(number,ax)
    def vec_b(self):
        # makes a vector of magnetic fields in the same ordering as
        # the_matrix class below
        the_vector=np.zeros((self.numsensors*3))
        for j in range(myarray.numsensors):
            r = myarray.sensors[j].pos
            b=np.array([bxtarget(r[0],r[1],r[2]),
                        bytarget(r[0],r[1],r[2]),
                        bztarget(r[0],r[1],r[2])])
            for k in range(3):
                the_vector[j*3+k]=b[k]
        return the_vector


# set up array of sensors
a_sensors=0.5
p0=np.array([-a_sensors/2,-a_sensors/2,-a_sensors/2])
p1=p0+np.array([a_sensors,0,0])
p2=p0+np.array([0,a_sensors,0])
p3=p0+np.array([0,0,a_sensors])
points=(p0,p1,p2,p3)

nsensors=int(options.nsensors)
myarray=sensorarray(nsensors,nsensors,nsensors,points)
print(myarray.sensors[0].pos)
print(myarray.numsensors)
print(myarray.sensors[myarray.numsensors-1].pos)
print(myarray.sensors[myarray.numsensors-2].pos)

print('the vector test')
print(len(myarray.vec_b()),myarray.vec_b())

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

if(options.traces):
    fig = plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    myset.draw_coils(ax)
    myarray.draw_sensors(ax)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()


from matplotlib import cm

class the_matrix:
    def __init__(self,myset,myarray):
        self.m=np.zeros((myset.numcoils,myarray.numsensors*3))
        #self.fill(myset,myarray)
        self.fillspeed(myset,myarray)
        self.condition = np.linalg.cond(self.m)

        # for some reason I chose to create the transpose of the usual
        # convention, when I first wrote the fill method
        self.capital_M=self.m.T # M=s*c=sensors*coils Matrix

        # Do the svd
        self.U,self.s,self.VT=np.linalg.svd(self.capital_M)

        print('s is',self.s)
        # s is just a list of the diagonal elements, rather than a matrix
        # You can make the matrix this way:
        self.S=np.zeros(self.capital_M.shape)
        self.S[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(self.s)
        # Or, I've seen people use "full_matrices=True" in the svd command

        # Start to calculate the inverse explicitly
        # list of reciprocals
        d=1./self.s
        self.D=np.zeros(self.capital_M.shape)
        # matrix of reciprocals
        self.D[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(d)

        # inverse of capital_M
        self.Minv=self.VT.T.dot(self.D.T).dot(self.U.T)
        #self.Minv=np.linalg.pinv(self.capital_M)
        
        # now gets to fixin'
        # remove just the last mode
        n_elements=myset.numcoils-1

        self.Dp=self.D[:,:n_elements]
        self.VTp=self.VT[:n_elements,:]
        self.Minvp=self.VTp.T.dot(self.Dp.T).dot(self.U.T)
        
    def fill(self,myset,myarray):
        for i in range(myset.numcoils):
            myset.set_independent_current(i,1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b = myset.b(r)
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
            myset.set_independent_current(i,0.0)

    def fillspeed(self,myset,myarray):
        myset.set_common_current(1.0)
        for i in range(myset.numcoils):
            print("Doing coil %d"%i)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                bx,by,bz=myset.coil[i].b_prime(r[0],r[1],r[2])
                b=[bx,by,bz]
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
        myset.zero_currents()
            
    def check_field_graphically(self,myset,myarray):
        # test each coil by graphing field at each sensor
        for i in range(myset.numcoils):
            fig = plt.figure()
            ax=fig.add_subplot(111,projection='3d')
            myset.draw_coil(i,ax)
            myset.coil[i].set_current(1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b=myset.b(r)
                bhat=b*5.e4
                points = []
                points.append(r)
                points.append(r+bhat)
                xs = ([p[0] for p in points])
                ys = ([p[1] for p in points])
                zs = ([p[2] for p in points])
                ax.plot(xs,ys,zs)
            myset.coil[i].set_current(0.0)
            ax.legend()
            plt.show()

    def show_matrices(self):
        fig1,ax1=plt.subplots()
        fig2,ax2=plt.subplots()
        fig3,ax3=plt.subplots()
        fig4,ax4=plt.subplots()
        fig5,ax5=plt.subplots()
        fig6,ax6=plt.subplots()
        fig7,ax7=plt.subplots()
        fig8,ax8=plt.subplots()
        fig9,ax9=plt.subplots()

        ax1.imshow(self.capital_M,cmap=cm.bwr)
        ax2.imshow(self.U,cmap=cm.bwr)
        ax3.imshow(self.S,cmap=cm.bwr)
        ax4.imshow(self.VT,cmap=cm.bwr)
        ax5.imshow(self.D,cmap=cm.bwr)
        ax6.imshow(self.Minv,cmap=cm.bwr)
        ax7.imshow(self.Dp,cmap=cm.bwr)
        ax8.imshow(self.VTp,cmap=cm.bwr)
        ax9.imshow(self.Minvp,cmap=cm.bwr)

        ax1.set_title('M')
        ax2.set_title('U')
        ax3.set_title('S')
        ax4.set_title('VT')
        ax5.set_title('D')
        ax6.set_title('Minv')
        ax7.set_title('Dp')
        ax8.set_title('VTp')
        ax9.set_title('Minvp')

        plt.show()

        
mymatrix=the_matrix(myset,myarray)

print('The condition number is %f'%mymatrix.condition)
if(options.matrices):
    mymatrix.show_matrices()

# Set up vector of desired fields

#print(len(myarray.vec_b()),myarray.vec_b())
#vec_i=mymatrix.Minvp.dot(myarray.vec_b())
vec_i=mymatrix.Minv.dot(myarray.vec_b())
#print(vec_i)

# Assign currents to coilcube

myset.set_currents(vec_i)

# Check the field at the center of the coilcube
r=np.array([0,0,0])
print(myset.b(r))
print(myset.b_prime(0,0,0))

from scipy.optimize import curve_fit

def fiteven(x,p0,p2,p4,p6):
    return p0+p2*x**2+p4*x**4+p6*x**6

def fitodd(x,p1,p3,p5,p7):
    return p1*x+p3*x**3+p5*x**5+p7*x**7

def fitgraph(xdata,ydata,ax):
    popt,pcov=curve_fit(fiteven,xdata[abs(xdata)<.5],ydata[abs(xdata)<.5])
    print(popt)
    ax.plot(points1d,fiteven(xdata,*popt),'r--',label='$p_0$=%2.1e,$p_2$=%2.1e,$p_4$=%2.1e,$p_6$=%2.1e'%tuple(popt))

# scans along each axis
points1d=np.mgrid[-1:1:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=myset.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=myset.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=myset.b_prime(0.,0.,points1d)

# target field
bx1d_target_xscan=bxtarget(points1d,0.,0.)*np.ones(np.shape(points1d))
bx1d_target_yscan=bxtarget(0.,points1d,0.)*np.ones(np.shape(points1d))
bx1d_target_zscan=bxtarget(0.,0.,points1d)*np.ones(np.shape(points1d))

by1d_target_xscan=bytarget(points1d,0.,0.)*np.ones(np.shape(points1d))
by1d_target_yscan=bytarget(0.,points1d,0.)*np.ones(np.shape(points1d))
by1d_target_zscan=bytarget(0.,0.,points1d)*np.ones(np.shape(points1d))

bz1d_target_xscan=bztarget(points1d,0.,0.)*np.ones(np.shape(points1d))
bz1d_target_yscan=bztarget(0.,points1d,0.)*np.ones(np.shape(points1d))
bz1d_target_zscan=bztarget(0.,0.,points1d)*np.ones(np.shape(points1d))

if(options.zoom):
    mask=(points1d>=-a_sensors/2)&(points1d<=a_sensors/2)
else:
    mask=np.full(np.shape(points1d),True)

if(options.axes):
    fig7,(ax71)=plt.subplots(nrows=1)
    fig8,(ax81)=plt.subplots(nrows=1)
    fig9,(ax91)=plt.subplots(nrows=1)
    
    ax71.plot(points1d[mask],bz1d_xscan[mask],label='$B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_target_xscan[mask],label='target $B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_yscan[mask],label='$B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_target_yscan[mask],label='target $B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_zscan[mask],label='$B_z(0,0,z)$')
    ax71.plot(points1d[mask],bz1d_target_zscan[mask],label='target $B_z(0,0,z)$')
    ax71.set_xlabel('x, y, or z (m)')
    from sympy import latex
    if(options.dipole):
        ax71.set_ylabel('$B_z=dipole$')
    else:
        ax71.set_ylabel('$B_z=\Pi_{z,%d,%d}=%s$'%(l,m,latex(sp.Piz)))
    if(not options.zoom):
        ax71.axvline(x=a/2,color='black',linestyle='--')
        ax71.axvline(x=-a/2,color='black',linestyle='--')
        ax71.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax71.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax81.plot(points1d[mask],by1d_xscan[mask],label='$B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_target_xscan[mask],label='target $B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_yscan[mask],label='$B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_target_yscan[mask],label='target $B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_zscan[mask],label='$B_y(0,0,z)$')
    ax81.plot(points1d[mask],by1d_target_zscan[mask],label='target $B_y(0,0,z)$')
    ax81.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax81.set_ylabel('$B_y=dipole$')
    else:
        ax81.set_ylabel('$B_y=\Pi_{y,%d,%d}=%s$'%(l,m,latex(sp.Piy)))
    if(not options.zoom):
        ax81.axvline(x=a/2,color='black',linestyle='--')
        ax81.axvline(x=-a/2,color='black',linestyle='--')
        ax81.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax81.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax91.plot(points1d[mask],bx1d_xscan[mask],label='$B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_target_xscan[mask],label='target $B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_yscan[mask],label='$B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_target_yscan[mask],label='target $B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_zscan[mask],label='$B_x(0,0,z)$')
    ax91.plot(points1d[mask],bx1d_target_zscan[mask],label='target $B_x(0,0,z)$')
    ax91.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax91.set_ylabel('$B_x=dipole$')
    else:
        ax91.set_ylabel('$B_x=\Pi_{x,%d,%d}=%s$'%(l,m,latex(sp.Pix)))
    if(not options.zoom):
        ax91.axvline(x=a/2,color='black',linestyle='--')
        ax91.axvline(x=-a/2,color='black',linestyle='--')
        ax91.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax91.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax71.axhline(y=0,color='black')
    ax81.axhline(y=0,color='black')
    ax91.axhline(y=0,color='black')
    
    ax71.legend()
    ax81.legend()
    ax91.legend()


if(options.residuals):

    ax101=plt.figure(101)
    plt.plot(points1d[mask],bz1d_xscan[mask]-bz1d_target_xscan[mask],label='residual $B_z(x,0,0)$')
    plt.plot(points1d[mask],bz1d_yscan[mask]-bz1d_target_yscan[mask],label='residual $B_z(0,y,0)$')
    plt.plot(points1d[mask],bz1d_zscan[mask]-bz1d_target_zscan[mask],label='residual $B_z(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_z$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax102=plt.figure(102)
    plt.plot(points1d[mask],by1d_xscan[mask]-by1d_target_xscan[mask],label='residual $B_y(x,0,0)$')
    plt.plot(points1d[mask],by1d_yscan[mask]-by1d_target_yscan[mask],label='residual $B_y(0,y,0)$')
    plt.plot(points1d[mask],by1d_zscan[mask]-by1d_target_zscan[mask],label='residual $B_y(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_y$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax103=plt.figure(103)
    plt.plot(points1d[mask],bx1d_xscan[mask]-bx1d_target_xscan[mask],label='residual $B_x(x,0,0)$')
    plt.plot(points1d[mask],bx1d_yscan[mask]-bx1d_target_yscan[mask],label='residual $B_x(0,y,0)$')
    plt.plot(points1d[mask],bx1d_zscan[mask]-bx1d_target_zscan[mask],label='residual $B_x(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_x$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

plt.show()

# studies over an ROI

#x,y,z=np.mgrid[-.25:.25:51j,-.25:.25:51j,-.25:.25:51j]
#x,y,z=np.mgrid[-.49:.49:99j,-.49:.49:99j,-.49:.49:99j]
x,y,z=np.mgrid[-.5:.5:101j,-.5:.5:101j,-.5:.5:101j]

if(options.incells):
    rcell=0.3 # m, cell radius
    hcell=0.1601 # m, cell height
    dcell=0.08 # m, bottom to top distance of cells
    mask=(abs(z)>=dcell/2)&(abs(z)<=dcell/2+hcell)&(x**2+y**2<rcell**2)
    mask_upper=(abs(z)>=dcell/2)&(abs(z)<=dcell/2+hcell)&(x**2+y**2<rcell**2)&(z>0)
    mask_lower=(abs(z)>=dcell/2)&(abs(z)<=dcell/2+hcell)&(x**2+y**2<rcell**2)&(z<0)
else:
    mask=np.full(np.shape(z),True)
    mask_upper=(z>0)
    mask_lower=(z<0)

# This is used to test the cell dimensions.

#fig=plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#scat=ax.scatter(x[mask_upper],y[mask_upper],z[mask_upper])
#plt.show()

bx_roi,by_roi,bz_roi=myset.b_prime(x,y,z)
bx_target=bxtarget(x,y,z)
by_target=bytarget(x,y,z)
bz_target=bztarget(x,y,z)
bx_residual=bx_roi-bx_target
by_residual=by_roi-by_target
bz_residual=bz_roi-bz_target

print(np.shape(bx_roi))

print('Statistics on the ROI')
print



bz_ave=np.average(bz_target)
print('The unmasked average Bz prior to correction is %e'%bz_ave)
bz_max=np.amax(bz_target)
bz_min=np.amin(bz_target)
bz_delta=bz_max-bz_min
print('The unmasked max/min/diff Bz are %e %e %e'%(bz_max,bz_min,bz_delta))
print('We normalize this to 3 nT max-min')
print

print('Both cells')
bz_mask_max=np.amax(bz_target[mask])
bz_mask_min=np.amin(bz_target[mask])
bz_mask_delta=bz_mask_max-bz_mask_min
print('The max/min/diff Bz masks are %e %e %e'%(bz_mask_max,bz_mask_min,bz_mask_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_mask_delta/bz_delta*3))
bz_std=np.std(bz_target[mask])
print('The masked standard deviation of Bz is %e'%bz_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_std/bz_delta*3))
print

bz_residual_max=np.amax(bz_residual[mask])
bz_residual_min=np.amin(bz_residual[mask])
bz_residual_delta=bz_residual_max-bz_residual_min
print('The max/min/diff Bz residuals are %e %e %e'%(bz_residual_max,bz_residual_min,bz_residual_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_residual_delta/bz_delta*3))
bz_residual_std=np.std(bz_residual[mask])
print('The standard deviation of Bz residuals is %e'%bz_residual_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_residual_std/bz_delta*3))
print

print('Upper cell')
bz_mask_max=np.amax(bz_target[mask_upper])
bz_mask_min=np.amin(bz_target[mask_upper])
bz_mask_delta=bz_mask_max-bz_mask_min
print('The max/min/diff Bz masks are %e %e %e'%(bz_mask_max,bz_mask_min,bz_mask_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_mask_delta/bz_delta*3))
bz_std=np.std(bz_target[mask_upper])
print('The masked standard deviation of Bz is %e'%bz_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_std/bz_delta*3))
print

bz_residual_max=np.amax(bz_residual[mask_upper])
bz_residual_min=np.amin(bz_residual[mask_upper])
bz_residual_delta=bz_residual_max-bz_residual_min
print('The max/min/diff Bz residuals are %e %e %e'%(bz_residual_max,bz_residual_min,bz_residual_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_residual_delta/bz_delta*3))
bz_residual_std=np.std(bz_residual[mask_upper])
print('The standard deviation of Bz residuals is %e'%bz_residual_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_residual_std/bz_delta*3))
print

print('Lower cell')
bz_mask_max=np.amax(bz_target[mask_lower])
bz_mask_min=np.amin(bz_target[mask_lower])
bz_mask_delta=bz_mask_max-bz_mask_min
print('The max/min/diff Bz masks are %e %e %e'%(bz_mask_max,bz_mask_min,bz_mask_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_mask_delta/bz_delta*3))
bz_std=np.std(bz_target[mask_lower])
print('The masked standard deviation of Bz is %e'%bz_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_std/bz_delta*3))
print

bz_residual_max=np.amax(bz_residual[mask_lower])
bz_residual_min=np.amin(bz_residual[mask_lower])
bz_residual_delta=bz_residual_max-bz_residual_min
print('The max/min/diff Bz residuals are %e %e %e'%(bz_residual_max,bz_residual_min,bz_residual_delta))
print('Normalizing to 3 nT gives a delta of %f nT'%(bz_residual_delta/bz_delta*3))
bz_residual_std=np.std(bz_residual[mask_lower])
print('The standard deviation of Bz residuals is %e'%bz_residual_std)
print('Normalizing to 3 nT gives a standard deviation of %f nT'%(bz_residual_std/bz_delta*3))
print


bt2_target=bx_target**2+by_target**2+bz_target**2
bt2_ave=np.average(bt2_target[mask])
print('The BT2 prior to correction is %e'%bt2_ave)
bt2_ave_norm=bt2_ave*3**2/bz_delta**2
print('Normalized is %f nT^2'%(bt2_ave_norm))

bt2_residual=bx_residual**2+by_residual**2+bz_residual**2
bt2_residual_ave=np.average(bt2_residual[mask])
print('The BT2 after correction is %e'%bt2_residual_ave)
bt2_residual_ave_norm=bt2_residual_ave*3**2/bz_delta**2
print('Normalized is %f nT^2'%(bt2_residual_ave_norm))
print

print('The normalized currents are:')
vec_i=vec_i*3e-9/bz_delta
print(vec_i)
print('The maximum current is %f A'%np.amax(vec_i))
print('The minimum current is %f A'%np.amin(vec_i))

n=14 # number of bits
#Igoal=round(np.amax(vec_i)*1000)
Imin=-0.025 # minimum current (A)
Imax=0.025 # maximum current (A)
deltaI=Imax-Imin # full range of current, which I'll distribute my bits across.

print("Distributing %d bits across %f A"%(n,deltaI))

def bits(I):
    return np.rint((I-Imin)/deltaI*(2**n-1))

def true_current(b):
    return deltaI/(2**n-1)*b+Imin

print('The bits are',bits(vec_i))

dig_i=true_current(bits(vec_i))

print('The digitized currents are',dig_i)

dig_i_reunnormalized=dig_i/3e-9*bz_delta

#trying to plot normalized and digitized currents vs number coils? 

#cn=list(range(0,50)) #number of coils
#plt.plot(cn,vec_i,label='normalized currents')
#plt.plot(cn,true_current(bits(vec_i)),label='digitized currents')
#plt.show()

myset.set_currents(dig_i_reunnormalized)

points1d=np.mgrid[-1:1:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=myset.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=myset.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=myset.b_prime(0.,0.,points1d)

if(options.zoom):
    mask=(points1d>=-a_sensors/2)&(points1d<=a_sensors/2)
else:
    mask=np.full(np.shape(points1d),True)

if(options.axes):
    fig7,(ax71)=plt.subplots(nrows=1)
    fig8,(ax81)=plt.subplots(nrows=1)
    fig9,(ax91)=plt.subplots(nrows=1)
    
    ax71.plot(points1d[mask],bz1d_xscan[mask],label='$B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_target_xscan[mask],label='target $B_z(x,0,0)$')
    ax71.plot(points1d[mask],bz1d_yscan[mask],label='$B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_target_yscan[mask],label='target $B_z(0,y,0)$')
    ax71.plot(points1d[mask],bz1d_zscan[mask],label='$B_z(0,0,z)$')
    ax71.plot(points1d[mask],bz1d_target_zscan[mask],label='target $B_z(0,0,z)$')
    ax71.set_xlabel('x, y, or z (m)')
    from sympy import latex
    if(options.dipole):
        ax71.set_ylabel('$B_z=dipole$')
    else:
        ax71.set_ylabel('$B_z=\Pi_{z,%d,%d}=%s$'%(l,m,latex(sp.Piz)))
    if(not options.zoom):
        ax71.axvline(x=a/2,color='black',linestyle='--')
        ax71.axvline(x=-a/2,color='black',linestyle='--')
        ax71.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax71.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax81.plot(points1d[mask],by1d_xscan[mask],label='$B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_target_xscan[mask],label='target $B_y(x,0,0)$')
    ax81.plot(points1d[mask],by1d_yscan[mask],label='$B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_target_yscan[mask],label='target $B_y(0,y,0)$')
    ax81.plot(points1d[mask],by1d_zscan[mask],label='$B_y(0,0,z)$')
    ax81.plot(points1d[mask],by1d_target_zscan[mask],label='target $B_y(0,0,z)$')
    ax81.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax81.set_ylabel('$B_y=dipole$')
    else:
        ax81.set_ylabel('$B_y=\Pi_{y,%d,%d}=%s$'%(l,m,latex(sp.Piy)))
    if(not options.zoom):
        ax81.axvline(x=a/2,color='black',linestyle='--')
        ax81.axvline(x=-a/2,color='black',linestyle='--')
        ax81.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax81.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax91.plot(points1d[mask],bx1d_xscan[mask],label='$B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_target_xscan[mask],label='target $B_x(x,0,0)$')
    ax91.plot(points1d[mask],bx1d_yscan[mask],label='$B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_target_yscan[mask],label='target $B_x(0,y,0)$')
    ax91.plot(points1d[mask],bx1d_zscan[mask],label='$B_x(0,0,z)$')
    ax91.plot(points1d[mask],bx1d_target_zscan[mask],label='target $B_x(0,0,z)$')
    ax91.set_xlabel('x, y, or z (m)')
    if(options.dipole):
        ax91.set_ylabel('$B_x=dipole$')
    else:
        ax91.set_ylabel('$B_x=\Pi_{x,%d,%d}=%s$'%(l,m,latex(sp.Pix)))
    if(not options.zoom):
        ax91.axvline(x=a/2,color='black',linestyle='--')
        ax91.axvline(x=-a/2,color='black',linestyle='--')
        ax91.axvline(x=a_sensors/2,color='red',linestyle='--')
        ax91.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax71.axhline(y=0,color='black')
    ax81.axhline(y=0,color='black')
    ax91.axhline(y=0,color='black')
    
    ax71.legend()
    ax81.legend()
    ax91.legend()


if(options.residuals):

    ax101=plt.figure(101)
    plt.plot(points1d[mask],bz1d_xscan[mask]-bz1d_target_xscan[mask],label='residual $B_z(x,0,0)$')
    plt.plot(points1d[mask],bz1d_yscan[mask]-bz1d_target_yscan[mask],label='residual $B_z(0,y,0)$')
    plt.plot(points1d[mask],bz1d_zscan[mask]-bz1d_target_zscan[mask],label='residual $B_z(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_z$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax102=plt.figure(102)
    plt.plot(points1d[mask],by1d_xscan[mask]-by1d_target_xscan[mask],label='residual $B_y(x,0,0)$')
    plt.plot(points1d[mask],by1d_yscan[mask]-by1d_target_yscan[mask],label='residual $B_y(0,y,0)$')
    plt.plot(points1d[mask],by1d_zscan[mask]-by1d_target_zscan[mask],label='residual $B_y(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_y$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

    ax103=plt.figure(103)
    plt.plot(points1d[mask],bx1d_xscan[mask]-bx1d_target_xscan[mask],label='residual $B_x(x,0,0)$')
    plt.plot(points1d[mask],bx1d_yscan[mask]-bx1d_target_yscan[mask],label='residual $B_x(0,y,0)$')
    plt.plot(points1d[mask],bx1d_zscan[mask]-bx1d_target_zscan[mask],label='residual $B_x(0,0,z)$')
    plt.xlabel('x, y, or z (m)')
    plt.ylabel('residual $B_x$ (true-target)')
    plt.legend()
    if(not options.zoom):
        plt.axvline(x=a/2,color='black',linestyle='--')
        plt.axvline(x=-a/2,color='black',linestyle='--')
        plt.axvline(x=a_sensors/2,color='red',linestyle='--')
        plt.axvline(x=-a_sensors/2,color='red',linestyle='--')

plt.show()
