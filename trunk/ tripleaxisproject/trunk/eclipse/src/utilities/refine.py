from __future__ import division
import numpy as N
import rescalculator.lattice_calculator as lattice_calculator
pi=N.pi

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
    return astar[0],alphastar[0]
    
    

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
    print 'alphastar',alphastar
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
    print 'Q',Q
    d=2*pi/Q
    print 'd',d
    return d
    
 
def calc_cosn(h,k,l,d,alphastar,astar,phi):
    """These H,K,L,should be in rhombohedral"""
    n=h**2+k**2+l**2
    r=h*k+k*l+l*h
    cosn=(n-r)*(1-cos(alphastar))*N.sin(phi)**2+\
        (n+2*r)*(1+2*cos(alphastar))*N.cos(phi)**2
    
    cosn=cosn*astar**2*d**2/4/pi**2/3

def pseudocubic2hex(h,k,l):
    g=N.array([[1,-1,0],[0,1,-1],[2,2,2]],'Float64')
    hh=[]
    kh=[]
    lh=[]
    for i in range(len(l)):
        vec=N.array([h[i],k[i],l[i]],'Float64')
        vech=N.dot(g,vec)
        hh.append(vech[0])
        kh.append(vech[1])
        lh.append(vech[2])
        
    return N.array(hh),N.array(kh),N.array(lh)
    
    
def hex2rhomb(h,k,l):
    g=N.array([[2,1,1],[-1,1,1],[-1,-2,1]],'Float64')/3
    hr=[]
    kr=[]
    lr=[]
    for i in range(len(l)):
        vec=N.array([h[i],k[i],l[i]],'Float64')
        vecr=N.dot(g,vec)
        hr.append(vecr[0])
        kr.append(vecr[1])
        lr.append(vecr[2])
        
    return N.array(hr),N.array(kr),N.array(lr)

    
if __name__=="__main__":
    Hpc=N.array([1,1,2,1],'float64')
    Kpc=N.array([1,1,-1,0],'float64')
    Lpc=N.array([1,-2,-1,0],'float64')
    Hh,Kh,Lh=pseudocubic2hex(Hpc,Kpc,Lpc)
    print 'hex'
    for i in range(len(Hh)):
        print Hh[i],Kh[i],Lh[i]
    Hr,Kr,Lr=hex2rhomb(Hh,Kh,Lh)
    print 'rhomb'
    for i in range(len(Hr)):
        print Hr[i],Kr[i],Lr[i]
    #Hh=N.array([0,1,1,2],'float64')
    #Kh=N.array([1,1,1,-1],'float64')
    #Lh=N.array([2,1,-2,-1],'float64')
    #d=calcd(Hh,Kh,Lh)
    #astar,alphastar=calcstar()
    #print astar,N.degrees(alphastar)
    
    
    
    