#!/usr/bin/python

# Fri May 24 10:11:45 CDT 2019 Jeff added this line.


from scipy.constants import mu_0, pi
import numpy as np

def b_segment(i,p0,p1,r):
    # p0 is one end (vector in m)
    # p1 is the other (m)
    # r is the position of interest (m)
    # i is the current (A)
    d0 = r - p0
    d1 = r - p1
    ell = p1 - p0
    lend0 = np.sqrt(d0.dot(d0))
    lend1 = np.sqrt(d1.dot(d1))
    lenell = np.sqrt(ell.dot(ell))
    costheta0 = np.inner(ell,d0)/lenell/lend0
    costheta1 = -np.inner(ell,d1)/lenell/lend1
    ellcrossd0 = np.cross(ell,d0)
    lenellcrossd0 = np.sqrt(ellcrossd0.dot(ellcrossd0))
    modsintheta0 = lenellcrossd0/lenell/lend0
    a = lend0 * modsintheta0
    if(lenellcrossd0>0):
        nhat = ellcrossd0/lenellcrossd0
    else:
        nhat = np.array([0,0,0])

    if(a>0):
        b_total=mu_0*i/4.0/pi/a*(costheta0+costheta1)*nhat
    else:
        b_total = np.array([0,0,0])

    return b_total

def b_loop(i,points,r):
    # i is the current (A)
    # points is a list of numpy 3-arrays defining the loop (m)
    # r is the position of interest (m)
    b_total = np.array([0,0,0])
    for j in range(len(points)):
        b_total = b_total + b_segment(i,points[j-1],points[j],r)
    return b_total

# Halliday & Resnick, 10th ed., question 29.13
p0 = np.array([0,0,0])
p1 = np.array([0.18,0,0])
r = np.array([0.09,0.131,0])
i = 0.0582
print(b_segment(i,p0,p1,r))

# Halliday & Resnick, 10th ed., question 29.17
p0 = np.array([0,0,0])
p1 = np.array([0.136,0,0])
r = np.array([0.136,0.251,0])
i = 0.693
print(b_segment(i,p0,p1,r))

# Halliday & Resnick, 10th ed., question 29.31
a = 0.047
i = 13.0
p0 = np.array([0,0,0])
p1 = np.array([2*a,0,0])
p2 = np.array([2*a,a,0])
p3 = np.array([a,a,0])
p4 = np.array([a,2*a,0])
p5 = np.array([0,2*a,0])
points = (p0,p1,p2,p3,p4,p5)
r = np.array([2*a,2*a,0])
print(b_loop(i,points,r))

# Halliday & Resnick, 10th ed., question 29.83
a = 0.08
i = 10.0
p0 = np.array([0,0,0])
p1 = np.array([a,0,0])
p2 = np.array([a,-a,0])
p3 = np.array([0,-a,0])
points = (p0,p1,p2,p3)
r = np.array([a/4,-a/4,0])
print(b_loop(i,points,r))

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
        self.numcoils = (xdim*ydim + xdim*zdim + ydim*zdim)*2
        self.numindep = self.numcoils - 1 # valid for xdim=ydim=zdim=1
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
        # wire the last two coils together
        # only works if xdim=ydim=zdim=1
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
            self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim].set_current(current)
    def draw_coil(self,number,ax):
        coil = self.coil(number)
        points = coil.corners + (coil.corners[0],)
        x = ([p[0] for p in points])
        y = ([p[1] for p in points])
        z = ([p[2] for p in points])
        ax.plot(x,y,z,label='coil')
    def draw_coils(self,ax):
        for number in range(self.numcoils):
            self.draw_coil(number,ax)
    def b(self,r):
        b_total = 0.0
        for number in range(self.numcoils):
            b_total = b_total + self.coil(number).b(r)
        return b_total

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

class coil:
    def __init__(self,corners,current):
        self.corners = corners
        self.current = current
    def set_current(self,current):
        self.current = current
    def b(self,r):
        return b_loop(self.current,self.corners,r)

# repeat of 29.83:
thiscoil = coil(points,i)
print(thiscoil.b(r))
thiscoil.set_current(i/2.0)
print(thiscoil.b(r))

# test of face class -- made a slight change after this test

# p0 is kind of like the origin; p1 and p3 are kind of like basis
# vectors, defining the sides of the rectangular face, whose corner is
# p0
points = (p0,p1,p3)
thisface = face(2,2,points)
print(thisface.coilnum)
print(thisface.coil[0].corners)
print(thisface.coil[1].corners)
print(thisface.coil[2].corners)
print(thisface.coil[3].corners)

# test of coilcube class

# again p0 is like the origin; p1, p2, and p3 define the x, y, and z
# (or whatever order) sides of the cube
print("Coilcube test")
a = 1.0
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
print('hello')
print(points)
mycube = coilcube(3,3,3,points)
#print(mycube.face[1].coil[2].corners)
#print(mycube.coil(6).corners)
#print(mycube.coil(6).current)
print(mycube.numcoils)

class sensor:
    def __init__(self,pos):
        self.pos = pos

class sensorarray:
    def __init__(self,xdim,ydim,zdim,corners):
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
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

# test of sensorarray class
a = 0.5
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
myarray = sensorarray(3,3,3,points)
print(myarray.sensors[0].pos)
print(myarray.numsensors)
print(myarray.sensors[myarray.numsensors-1].pos)
print(myarray.sensors[myarray.numsensors-2].pos)

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
mycube.draw_coils(ax)
myarray.draw_sensors(ax)
ax.legend()
plt.show()


print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(1.0)
print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(0.0)


from math import log10
from matplotlib import cm
from decimal import Decimal

class the_matrix:
    def __init__(self,mycube,myarray):
        self.m = np.zeros((mycube.numindep,myarray.numsensors*3))
        self.fill(mycube,myarray)
        self.minv=np.linalg.pinv(self.m)
        self.condition = np.linalg.cond(self.m)
        self.capital_M=self.m.T #(Here, M=s*c=sensors*coils Matrix  )
        self.capital_M1=np.linalg.pinv(self.capital_M)

        self.U, self.V, self.Wt = np.linalg.svd(self.capital_M, full_matrices=True)
        self.W, self.Vinv, self.Ut = np.linalg.svd(self.capital_M1, full_matrices=True)
        self.Vmat=np.zeros([len(np.arange(self.m.shape[1])),len(np.arange(self.m.shape[0]))]) 
        for d in range (0,len(np.arange(self.m.shape[0])),1):
	    self.Vmat[d][d]=self.V[d]
        self.Vmat_T=self.Vmat.T
        self.L_V=[]
        for i in range (0,len(np.arange(self.m.shape[0])),1):
	    self.L_V.append(round(log10(self.V[i]),1))
        print "The Diagonal Matrix in log form is : ", self.L_V
    def fill(self,mycube,myarray):
        # fill m, units of nT/A
        for i in range(mycube.numindep):
            mycube.set_independent_current(i,1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b = mycube.b(r)
                for k in range(3):
                    self.m[i,j*3+k]=b[k]*1.e9 # convert to nT here
            mycube.set_independent_current(i,0.0)
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
    def show_matrix(self):
        plt.imshow(self.m,interpolation='none')
        plt.colorbar()
        plt.show()
    def show_inverse(self):
        plt.imshow(self.minv,interpolation='none')
        plt.colorbar()
        plt.show()
    def show_matrices(self):
        #mp=1000000
        #mp4=1000
        #mp6=100
        #mp4s='%.1E' % Decimal(1000)
        #mp6s='%.1E' % Decimal(100)
        #mps='%.1E' % Decimal(1000000)
        #column_labels = ['1x', '1y', '1z','2x', '2y', '2z','3x', '3y', '3z','4x', '4y', '4z','6x', '6y', '6z', '8x', '8y', '8z']
        #row_labels_all = ['X-(1)', 'X+(2)', 'Y-(3)', 'Y+(4)', 'Z-(5)', 'Z+(6)']
        #row_labels=[]
        #for d in range(len(np.arange(self.m.shape[0]))):
	#    row_labels.append(row_labels_all[d])

        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        fig4, ax4 = plt.subplots()
        fig6, ax6 = plt.subplots()


        ax1.imshow(self.m, interpolation='none', cmap=cm.bwr, aspect='auto' )
        ax2.imshow(self.capital_M1, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
        ax3.imshow(self.Vmat, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
        ax4.imshow(self.U, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
        ax6.imshow(self.Wt, interpolation='nearest', cmap=cm.bwr, aspect='auto' )

        #ax1.set_xticklabels(column_labels)
        #ax2.set_xticklabels(column_labels)

        ax1.set_xticks(np.arange(self.m.shape[1]))
        ax1.set_yticks(np.arange(self.m.shape[0]))
        ax1.set_xlabel('Fluxgate positions')
        ax1.set_ylabel('Coils')
        #ax1.set_yticklabels(row_labels)
        ax1.set_title('Matrix M* (nT/A) ('+str(len(np.arange(self.m.shape[0])))+'coils * '+str(len(np.arange(self.m.shape[1])))+'sensors)')

        ax2.set_xticks(np.arange(self.capital_M1.shape[1]))
        ax2.set_yticks(np.arange(self.capital_M1.shape[0]))
        ax2.set_xlabel('Fluxgate positions')
        ax2.set_ylabel('Coils')
        #ax2.set_yticklabels(row_labels)
        #ax2.set_title('Pseudoinverse of M (*'+str(mps)+') (A/nT) ('+str(len(np.arange(self.m.shape[0])))+'coils * '+str(len(np.arange(self.m.shape[1])))+'sensors)')

        ax3.set_xticks(np.arange(self.Vmat.shape[1]))
        ax3.set_yticks(np.arange(self.Vmat.shape[0]))
        ax3.set_title('V-Sqrt of eigenvalues of M*M & MM* ('+str(len(np.arange(self.m.shape[1])))+'* '+str(len(np.arange(self.m.shape[0])))+')')

        ax4.set_xticks(np.arange(self.U.shape[1]))
        ax4.set_yticks(np.arange(self.U.shape[0]))
        #ax4.set_title('U-Orthonormal eigenvectors(*'+str(mp4s)+') of MM* ('+str(len(np.arange(self.m.shape[1])))+'*'+str(len(np.arange(self.m.shape[1])))+')')
        
        ax6.set_xticks(np.arange(self.Wt.shape[1]))
        ax6.set_yticks(np.arange(self.Wt.shape[0]))
        #ax6.set_title('W*-Orthonormal eigenvectors(*'+str(mp6s)+') of M*M ('+str(len(np.arange(self.m.shape[0])))+'*'+str(len(np.arange(self.m.shape[0])))+')')

        #for c in range (0,len(np.arange(self.m.shape[0]))):
	#    for s in range (0,len(np.arange(self.m.shape[1]))):
	#	ax1.text(s, c, int(self.m[c][s]), va='center', ha='center', rotation=90)
	#	ax2.text(s, c, int(self.capital_M1[c][s]*mp), va='center', ha='center', rotation=90)
	#	ax3.text(c, s, int(self.Vmat[s][c]), va='center', ha='center')

        #for c in range (0,len(np.arange(self.U.shape[0]))):
	#    for s in range (0,len(np.arange(self.U.shape[1]))):
	#	ax4.text(s, c, int(self.U[c][s]*mp4), va='center', ha='center')

        #for c in range (0,len(np.arange(self.Wt.shape[0]))):
	#    for s in range (0,len(np.arange(self.Wt.shape[1]))):
	#	ax6.text(s, c, int(self.Wt[c][s]*mp6), va='center', ha='center')

        plt.show()

        
mymatrix = the_matrix(mycube,myarray)

print(mymatrix.condition)
mymatrix.show_matrices()

