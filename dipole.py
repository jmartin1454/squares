from scipy.constants import mu_0, pi
from numpy import sqrt

class dipole:
    def __init__(self,xd,yd,zd,mx,my,mz):
        self.xd=xd
        self.yd=yd
        self.zd=zd
        self.mx=mx
        self.my=my
        self.mz=mz
    def b(self,x,y,z):
        xprime=x-self.xd
        yprime=y-self.yd
        zprime=z-self.zd
        rprime=sqrt(xprime**2+yprime**2+zprime**2)
        mdotrprime=self.mx*xprime+self.my*yprime+self.mz*zprime
        bx=mu_0/(4*pi)*(3*xprime*mdotrprime/rprime**5-self.mx/rprime**3)
        by=mu_0/(4*pi)*(3*yprime*mdotrprime/rprime**5-self.my/rprime**3)
        bz=mu_0/(4*pi)*(3*zprime*mdotrprime/rprime**5-self.mz/rprime**3)
        return bx,by,bz
    def bx(self,x,y,z):
        xprime=x-self.xd
        yprime=y-self.yd
        zprime=z-self.zd
        rprime=sqrt(xprime**2+yprime**2+zprime**2)
        mdotrprime=self.mx*xprime+self.my*yprime+self.mz*zprime
        bx=mu_0/(4*pi)*(3*xprime*mdotrprime/rprime**5-self.mx/rprime**3)
        return bx
    def by(self,x,y,z):
        xprime=x-self.xd
        yprime=y-self.yd
        zprime=z-self.zd
        rprime=sqrt(xprime**2+yprime**2+zprime**2)
        mdotrprime=self.mx*xprime+self.my*yprime+self.mz*zprime
        by=mu_0/(4*pi)*(3*yprime*mdotrprime/rprime**5-self.my/rprime**3)
        return by
    def bz(self,x,y,z):
        xprime=x-self.xd
        yprime=y-self.yd
        zprime=z-self.zd
        rprime=sqrt(xprime**2+yprime**2+zprime**2)
        mdotrprime=self.mx*xprime+self.my*yprime+self.mz*zprime
        bz=mu_0/(4*pi)*(3*zprime*mdotrprime/rprime**5-self.mz/rprime**3)
        return bz
