
import numpy as N
import rescalculator.lattice_calculator as lattice_calculator
pi=N.pi
#from spinwaves.utilities.mpfit.mpfit import mpfit
import sys,os,copy
import pylab
from openopt import NLSP
from numpy import cos, sin

A=5.581
B=A
C=13.8757


def calcstar():
    ar=N.sqrt(3*A**2+C**2)/3
    alphar=2*N.arcsin(3/2/N.sqrt(3+(C/A)**2))
    a=N.array([ar],'Float64')
    b=N.array([ar],'Float64')
    c=N.array([ar],'Float64')  
    alpha=N.array([alphar],'Float64')
    beta=N.array([alphar],'Float64')
    gamma=N.array([alphar],'Float64')
    orient1=N.array([[1,0,0]],'Float64')
    orient2=N.array([[0,0,1]],'Float64')
    orientation=lattice_calculator.Orientation(orient1,orient2)
    lattice = lattice_calculator.Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                 orientation=orientation)
    astar=lattice.astar
    alphastar=lattice.alphastar
    return astar[0],alphastar[0],lattice
    
    

def calcd(H,K,L):
    a=N.array([A],'Float64')
    b=N.array([B],'Float64')
    c=N.array([C],'Float64')
    
    
    alpha=N.radians(N.array([90],'Float64'))
    beta=N.radians(N.array([90],'Float64'))
    gamma=N.radians(N.array([120],'Float64'))
    orient1=N.array([[1,0,0]],'Float64')
    orient2=N.array([[0,0,1]],'Float64')
    orientation=lattice_calculator.Orientation(orient1,orient2)
    lattice = lattice_calculator.Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                                 orientation=orientation)
    newinput=lattice_calculator.CleanArgs(a=lattice.a,b=lattice.b,c=lattice.c,alpha=lattice.alpha,beta=lattice.beta,\
                             gamma=lattice.gamma,orient1=lattice._orient1,orient2=lattice._orient2,H=H,K=K,L=L)
    orientation=lattice_calculator.Orientation(newinput['orient1'],newinput['orient2'])
    lattice.__init__(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                        beta=newinput['beta'],gamma=newinput['gamma'],orientation=orientation\
                        )
    alphastar=lattice.alphastar
    print('alphastar',alphastar)
    EXP={}
    EXP['ana']={}
    EXP['ana']['tau']='pg(002)'
    EXP['mono']={}
    EXP['mono']['tau']='pg(002)';
    EXP['ana']['mosaic']=30
    EXP['mono']['mosaic']=30
    EXP['sample']={}
    EXP['sample']['mosaic']=10
    EXP['sample']['vmosaic']=10
    EXP['hcol']=N.array([40, 10, 20, 80],'Float64')
    EXP['vcol']=N.array([120, 120, 120, 120],'Float64')
    EXP['infix']=-1 #positive for fixed incident energy
    EXP['efixed']=14.7
    EXP['method']=0
    setup=[EXP]
    qx,qy,qz,Q=lattice.R2S(H,K,L)
    print('Q',Q)
    d=2*pi/Q
    print('d',d)
    return d

def setup(Hpc,Kpc,Lpc):
    Hh,Kh,Lh=pseudocubic2hex(Hpc,Kpc,Lpc)
    if 0:
        print('hex')
        for i in range(len(Hh)):
            print(Hh[i],Kh[i],Lh[i])
    Hr,Kr,Lr=hex2rhomb(Hh,Kh,Lh)    
    if 0:
        print('rhomb')
        for i in range(len(Hr)):
            print(Hr[i],Kr[i],Lr[i])
    d=calcd(Hh,Kh,Lh)
    astar,alphastar,lattice=calcstar()
    #print astar,N.degrees(alphastar)
    return Hr,Kr,Lr,d,astar,alphastar,lattice,Hh,Kh,Lh
    
def mgnfacFe3psquared(x):
    y=( 0.3972*N.exp(-13.244*x**2)+0.6295*N.exp(-4.903*x**2)\
   -0.0314*N.exp(+.35*x**2)-0.0044)**2
    return y

def calc_struct(p,h,k,l,d,q,alphastar,astar,lattice,hh,kh,lh):
    """h,k,l are rhombohedral, hh,kh,lh are hexagonal"""
    #Fe1h=N.array([0.0000,  0.0000,  0.2200],'Float64')
    #Fe2h=N.array([0.0000,  0.0000,  0.7200],'Float64')
    Fe1r=N.array([0.2200,  0.2200,  0.2200],'Float64')
    Fe2r=N.array([0.7200,  0.7200,  0.7200],'Float64')
    scale,theta1,phi1,theta2,phi2=p
    #cos2n1=calc_cos2n(phi1,h,k,l,d,alphastar,astar)
    #cos2n2=calc_cos2n(phi2,h,k,l,d,alphastar,astar)
    #s1=(1-cos2n1)
    #s2=(1-cos2n2)
    s1x=N.cos(theta1)*N.sin(phi1)
    s1y=N.sin(theta1)*N.sin(phi1)
    s1z=N.cos(phi1)
    s1=N.array([s1x,s1y,s1z],'Float64')
    s2x=N.cos(theta2)*N.sin(phi2)
    s2y=N.sin(theta2)*N.sin(phi2)
    s2z=N.cos(phi2)
    s2=N.array([s2x,s2y,s2z],'Float64')
    flist=[]
    #find fperpendicular from F-F.q*q hat
    for i in range(len(h)):
        #print 'i',i
        f=0
        vec=[hh[i],kh[i],lh[i]]
        dotp1=N.dot(vec,Fe1r)*2*pi
        dotp2=N.dot(vec,Fe2r)*2*pi
        sp1=s1-N.dot(s1,vec)*N.array(vec)/q[i]
        sp2=s2-N.dot(s2,vec)*N.array(vec)/q[i]
        f=sp1*N.exp(-1.0j*dotp1)+sp2*N.exp(-1.0j*dotp2)        
        f2=N.abs(N.dot(f,N.conjugate(f)))
        #print 'f1',f,f2
        
        vec=[kh[i],hh[i],lh[i]]
        dotp1=N.dot(vec,Fe1r)*2*pi
        dotp2=N.dot(vec,Fe2r)*2*pi
        f=sp1*N.exp(-1.0j*dotp1)+sp2*N.exp(-1.0j*dotp2)    
        f2=f2+N.abs(N.dot(f,N.conjugate(f)))
        #print 'f2',f,f2
        
        vec=[kh[i],-hh[i],lh[i]]
        dotp1=N.dot(vec,Fe1r)*2*pi
        dotp2=N.dot(vec,Fe2r)*2*pi
        f=sp1*N.exp(-1.0j*dotp1)+sp2*N.exp(-1.0j*dotp2)    
        f2=f2+N.abs(N.dot(f,N.conjugate(f)))
        #print 'f3',f,f2
        
        #f2=f2/3
        #print 'f2inal',f2
        flist.append(f2)
    flist=N.array(flist)*scale**2
    
    ff=mgnfacFe3psquared(q/4/pi)
    #print 'ff',ff
    Bf=.005
    flist=flist*ff*N.exp(-Bf*q**2/4/pi**2)
        
        
    return flist



def calcnew(p,v):
    th,a=p
    rot=a*N.array([[cos(th),-sin(th)],[sin(th),cos(th)]])
    res=N.dot(rot,v)
    return res

def firsteqn(p,v1,v2):
    res1=calcnew(p,v1)
    res2=calcnew(p,v2)
    outputvec=[res1[0]-2,
               res1[1]-0,
               res2[0]-2,
               res2[1]-2]
    return outputvec
    
    

if __name__=="__main__":
    th=N.radians(40.0)
    x0 = [th,1.0]
    v1=N.array([1.9564,-.00414])
    v2=N.array([1.93928,1.97075])
    p0 = NLSP(firsteqn, x0, args=(v1,v2))
    p0.ftol=1e-8
    p0.contol=1e-8
    r0 = p0.solve('nlp:ralg')
    print(N.degrees(r0.xf[0]))
    print(r0.xf[1])
    print(firsteqn(r0.xf,v1,v2))
    print('testing')
    print(calcnew(r0.xf,v1))
    print(calcnew(r0.xf,v2))
    print('testing')
    #Hr,Kr,Lr,d,astar,alphastar,lattice,Hh,Kh,Lh=setup(Hpc,Kpc,Lpc)
    
