#!/usr/bin/python3

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
# Define coil dimensions
import numpy as np

# Define coil dimensions and positions function
def calculate_coil_dimensions_and_positions(b,c,myset):
    # Calculate d
    d = (b - 3 * c) / 4
      
    
    # Calculate coil positions
    centre_list = []
    for i in range(3):
        x = -c - d + i * (c + d)
        centre_list.append(x)

    j=0
        
    for z in [-b/2,b/2]:
        for y in centre_list:
            for x in centre_list:
                point1=(x-c/2,y-c/2,z)
                point2=(x+c/2,y-c/2,z)
                point3=(x+c/2,y+c/2,z)
                point4=(x-c/2,y+c/2,z)
                myset.add_coil(np.array([point1, point2, point3, point4]))
                myset.set_current_in_coil(j,0.040)
                znew=z*(b/2+0.02)/(b/2) # calculate at 2 cm from centre
                print(znew)
                bx,by,bz=myset.b_prime(x,y,znew)
                print("Coil %d, Bx %e, By %e, Bz %e"%(j,bx,by,bz))
                myset.set_current_in_coil(j,0)
                j+=1

    for y in [-b/2,b/2]:
        for z in centre_list:
            for x in centre_list:
                point1=(x-c/2,y,z-c/2)
                point2=(x+c/2,y,z-c/2)
                point3=(x+c/2,y,z+c/2)
                point4=(x-c/2,y,z+c/2)
                myset.add_coil(np.array([point1, point2, point3, point4]))
                myset.set_current_in_coil(j,0.040)
                bx,by,bz=myset.b_prime(x,y,z)
                print("Coil %d, Bx %e, By %e, Bz %e"%(j,bx,by,bz))
                myset.set_current_in_coil(j,0)
                j+=1

    for x in [-b/2,b/2]:
        for z in centre_list:
            for y in centre_list:
                point1=(x,y-c/2,z-c/2)
                point2=(x,y+c/2,z-c/2)
                point3=(x,y+c/2,z+c/2)
                point4=(x,y-c/2,z+c/2)
                if not ((y==0)&((z==0)|(z==-c-d))):
                    myset.add_coil(np.array([point1, point2, point3, point4]))
                    myset.set_current_in_coil(j,0.040)
                    bx,by,bz=myset.b_prime(x,y,z)
                    print("Coil %d, Bx %e, By %e, Bz %e"%(j,bx,by,bz))
                    myset.set_current_in_coil(j,0)
                    j+=1


        
# Define box length and coil length

a=box_length =.26  # m
coil_length =.06  # m

# Calculate coil dimensions and positions
calculate_coil_dimensions_and_positions(box_length,coil_length,myset)


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
a_sensors=box_length/2
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
# the field at the center of the coilcube
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
points1d=np.mgrid[-a:a:101j]
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

print(vec_i)
max_unnormalized_current=np.max(np.abs(vec_i)) # arb. units
max_normalized_current=0.02 # Amperes
calibration_factor=max_normalized_current/max_unnormalized_current
calibrated_vec_i=vec_i*calibration_factor # Amperes
print(calibrated_vec_i)

# Now let's check what the field should be after setting these currents
myset.set_currents(calibrated_vec_i)

# the field at the centre of the coilcube
r=np.array([0,0,0])
print('Check field at centre of the coilcube')
print(myset.b(r))
print(myset.b_prime(0,0,0))


# scans along each axis
points1d=np.mgrid[-a:a:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=myset.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=myset.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=myset.b_prime(0.,0.,points1d)

# target field
bx1d_target_xscan=bxtarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
bx1d_target_yscan=bxtarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
bx1d_target_zscan=bxtarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor

by1d_target_xscan=bytarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
by1d_target_yscan=bytarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
by1d_target_zscan=bytarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor

bz1d_target_xscan=bztarget(points1d,0.,0.)*np.ones(np.shape(points1d))*calibration_factor
bz1d_target_yscan=bztarget(0.,points1d,0.)*np.ones(np.shape(points1d))*calibration_factor
bz1d_target_zscan=bztarget(0.,0.,points1d)*np.ones(np.shape(points1d))*calibration_factor

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

plt.show()
