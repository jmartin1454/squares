#!/usr/bin/python

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

parser.add_option("-c", "--ncoils", dest="ncoils", default=3,
                  help="nc where total coils is c=nc*nc*6")

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

d=dipole(1.2,0,0,0,0,1)

(options,args)=parser.parse_args()



class coilcube:
    def __init__(self,xdim,ydim,zdim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.zdim = zdim
        self.corners = corners
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        self.face = []
        thesecorners=(corners[0],corners[1],corners[2])
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=thesecorners+z
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=(corners[0],corners[1],corners[3])
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=thesecorners+y
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=(corners[0],corners[2],corners[3])
        self.face.append(face(ydim,zdim,thesecorners))
        thesecorners=thesecorners+x
        self.face.append(face(ydim,zdim,thesecorners))
        self.numcoils=(xdim*ydim+xdim*zdim+ydim*zdim)*2

    def coil(self,number):
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        if(number<xdim*ydim):
            return self.face[0].coil[number]
        elif(number<xdim*ydim*2):
            return self.face[1].coil[number-xdim*ydim]
        elif(number<xdim*ydim*2+xdim*zdim):
            return self.face[2].coil[number-xdim*ydim*2]
        elif(number<xdim*ydim*2+xdim*zdim*2):
            return self.face[3].coil[number-xdim*ydim*2-xdim*zdim]
        elif(number<xdim*ydim*2+xdim*zdim*2+ydim*zdim):
            return self.face[4].coil[number-xdim*ydim*2-xdim*zdim*2]
        else:
            return self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim]

    def set_independent_current(self,number,current):
        # # wire the last two coils together
        # # only works if xdim=ydim=zdim=1
        # made it back to not wired together Jeff
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        if(number<xdim*ydim):
            self.face[0].coil[number].set_current(current)
        elif(number<xdim*ydim*2):
            self.face[1].coil[number-xdim*ydim].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim):
            self.face[2].coil[number-xdim*ydim*2].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim*2):
            self.face[3].coil[number-xdim*ydim*2-xdim*zdim].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim*2+ydim*zdim):
            self.face[4].coil[number-xdim*ydim*2-xdim*zdim*2].set_current(current)
        else:
            self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim].set_current(current)

    def set_currents(self,vec_i):
        # set the currents to the vector given as the argument
        for i in range(self.numcoils):
            self.set_independent_current(i,vec_i[i])
            
    def zero_currents(self):
        # set all currents to zero
        for i in range(self.numcoils):
            self.set_independent_current(i,0.0)
            
    def set_common_current(self,curr):
        # set all currents to curr
        for i in range(self.numcoils):
            self.set_independent_current(i,curr)

    def draw_coil(self,number,ax):
        coil = self.coil(number)
        points = coil.points + (coil.points[0],)
        x = ([p[0] for p in points])
        y = ([p[1] for p in points])
        z = ([p[2] for p in points])
        ax.plot(x,y,z,label='coil')
        #ax.plot(x,y,z,label='coil')

    def draw_coils(self,ax):
        for number in range(self.numcoils):
            self.draw_coil(number,ax)

    def b(self,r):
        b_total = 0.0
        for number in range(self.numcoils):
            b_total = b_total + self.coil(number).b(r)
        return b_total

    def b_prime(self,x,y,z):
        b_total_x=0.*x
        b_total_y=0.*y
        b_total_z=0.*z
        for coilnum in range(self.numcoils):
            b_coil_x,b_coil_y,b_coil_z=self.coil(coilnum).b_prime(x,y,z)
            b_total_x=b_total_x+b_coil_x
            b_total_y=b_total_y+b_coil_y
            b_total_z=b_total_z+b_coil_z
        return b_total_x,b_total_y,b_total_z

class face:
    def __init__(self,xdim,ydim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.corners = corners
        x = corners[1]-corners[0]
        xstep = x/xdim
        y = corners[2]-corners[0]
        ystep = y/ydim
        coilnum = 0
        self.coil = []
        for i in range(xdim):
            for j in range(ydim):
                p0 = corners[0]+xstep*i+ystep*j
                p1 = p0+xstep
                p2 = p1+ystep
                p3 = p2-xstep
                points = (p0,p1,p2,p3)
                self.coil.append(coil(points,0.0))
                coilnum = coilnum + 1
        self.coilnum = coilnum


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

a=1.0
p0=np.array([-a/2,-a/2,-a/2])
p1=p0+np.array([a,0,0])
p2=p0+np.array([0,a,0])
p3=p0+np.array([0,0,a])
points=(p0,p1,p2,p3)

ncoils=int(options.ncoils)
mycube=coilcube(ncoils,ncoils,ncoils,points)

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
    ax = fig.gca(projection='3d')
    mycube.draw_coils(ax)
    myarray.draw_sensors(ax)
    #ax.legend()
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()


# print(mycube.b(myarray.sensors[0].pos))
# mycube.coil(0).set_current(1.0)
# print(mycube.b(myarray.sensors[0].pos))
# mycube.coil(0).set_current(0.0)

from matplotlib import cm

class the_matrix:
    def __init__(self,mycube,myarray):
        self.m=np.zeros((mycube.numcoils,myarray.numsensors*3))
        #self.fill(mycube,myarray)
        self.fillspeed(mycube,myarray)
        self.condition = np.linalg.cond(self.m)

        # for some reason I chose to create the transpose of the usual
        # convention, when I first wrote the fill method
        self.capital_M=self.m.T # M=s*c=sensors*coils Matrix

        # Do the svd
        self.U,self.s,self.VT=np.linalg.svd(self.capital_M)

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
        n_elements=mycube.numcoils-1

        self.Dp=self.D[:,:n_elements]
        self.VTp=self.VT[:n_elements,:]
        self.Minvp=self.VTp.T.dot(self.Dp.T).dot(self.U.T)
        
    def fill(self,mycube,myarray):
        for i in range(mycube.numcoils):
            mycube.set_independent_current(i,1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b = mycube.b(r)
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
            mycube.set_independent_current(i,0.0)

    def fillspeed(self,mycube,myarray):
        mycube.set_common_current(1.0)
        for i in range(mycube.numcoils):
            print("Doing coil %d"%i)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                bx,by,bz=mycube.coil(i).b_prime(r[0],r[1],r[2])
                b=[bx,by,bz]
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
        mycube.zero_currents()
            
    def check_field_graphically(self,mycube,myarray):
        # test each coil by graphing field at each sensor
        for i in range(mycube.numcoils):
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            mycube.draw_coil(i,ax)
            mycube.coil(i).set_current(1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b=mycube.b(r)
                bhat=b*5.e4
                points = []
                points.append(r)
                points.append(r+bhat)
                xs = ([p[0] for p in points])
                ys = ([p[1] for p in points])
                zs = ([p[2] for p in points])
                ax.plot(xs,ys,zs)
            mycube.coil(i).set_current(0.0)
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

        plt.show()

        
mymatrix = the_matrix(mycube,myarray)

print(mymatrix.condition)
if(options.matrices):
    mymatrix.show_matrices()

# Set up vector of desired fields

print(len(myarray.vec_b()),myarray.vec_b())
vec_i=mymatrix.Minvp.dot(myarray.vec_b())
print(vec_i)

# Assign currents to coilcube

mycube.set_currents(vec_i)

# Check the field at the center of the coilcube
r=np.array([0,0,0])
print(mycube.b(r))
print(mycube.b_prime(0,0,0))

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
bx1d_xscan,by1d_xscan,bz1d_xscan=mycube.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=mycube.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=mycube.b_prime(0.,0.,points1d)

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
x,y,z=np.mgrid[-.3:.3:61j,-.3:.3:61j,-.3:.3:61j]

if(options.incells):
    rcell=0.18 # m, cell radius
    hcell=0.15 # m, cell height
    dcell=0.10 # m, bottom to top distance of cells
    mask=(abs(z)>=dcell/2)&(abs(z)<=dcell/2+hcell)&(x**2+y**2<rcell**2)
else:
    mask=np.full(np.shape(z),True)

# This is used to test the cell dimensions.

#fig=plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#scat=ax.scatter(x[mask],y[mask],z[mask])
#plt.show()

bx_roi,by_roi,bz_roi=mycube.b_prime(x,y,z)
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
print('We normalize this to 0.6 nT max-min')
print

bz_std=np.std(bz_target[mask])
print('The masked standard deviation of Bz is %e'%bz_std)
print('Normalizing to 0.6 nT gives a standard deviation of %f nT'%(bz_std/bz_delta*0.6))
print

bz_residual_max=np.amax(bz_residual[mask])
bz_residual_min=np.amin(bz_residual[mask])
bz_residual_delta=bz_residual_max-bz_residual_min
print('The max/min/diff Bz residuals are %e %e %e'%(bz_residual_max,bz_residual_min,bz_residual_delta))
print('Normalizing to 0.6 nT gives a delta of %f nT'%(bz_residual_delta/bz_delta*0.6))
bz_residual_std=np.std(bz_residual[mask])
print('The standard deviation of Bz residuals is %e'%bz_residual_std)
print('Normalizing to 0.6 nT gives a standard deviation of %f nT'%(bz_residual_std/bz_delta*0.6))
print

bt2_target=bx_target**2+by_target**2+bz_target**2
bt2_ave=np.average(bt2_target[mask])
print('The BT2 prior to correction is %e'%bt2_ave)
bt2_ave_norm=bt2_ave*0.6**2/bz_delta**2
print('Normalized is %f nT^2'%(bt2_ave_norm))

bt2_residual=bx_residual**2+by_residual**2+bz_residual**2
bt2_residual_ave=np.average(bt2_residual[mask])
print('The BT2 after correction is %e'%bt2_residual_ave)
bt2_residual_ave_norm=bt2_residual_ave*0.6**2/bz_delta**2
print('Normalized is %f nT^2'%(bt2_residual_ave_norm))
