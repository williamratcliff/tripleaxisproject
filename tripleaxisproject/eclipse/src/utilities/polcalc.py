
import numpy as N
import rescalculator.lattice_calculator as lattice_calculator
pi=N.pi
#from spinwaves.utilities.mpfit.mpfit import mpfit
import sys,os,copy
import pylab
#from utilities.anneal import anneal

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
    
    

def calcQ(H,K,L):
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
    return Q
    
 
def calc_cos2n(phi,h,k,l,d,alphastar,astar):
    """These H,K,L,should be in rhombohedral"""
    n=h**2+k**2+l**2
    r=h*k+k*l+l*h
    cosn=(n-r)*(1-N.cos(alphastar))*N.sin(phi)**2+\
        (n+2*r)*(1+2*N.cos(alphastar))*N.cos(phi)**2
    
    cosn=cosn*astar**2*d**2/4/pi**2/3
    return cosn

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
    Q=calcQ(Hh,Kh,Lh)
    astar,alphastar,lattice=calcstar()
    #print astar,N.degrees(alphastar)
    return Q

def mgnfacFe3psquared(x):
    y=( 0.3972*N.exp(-13.244*x**2)+0.6295*N.exp(-4.903*x**2)\
   -0.0314*N.exp(+.35*x**2)-0.0044)**2
    return y

def genvec(vec):
    vmag=N.sqrt(vec[0]**2+vec[1]**2+vec[2]**2); vec=vec/vmag
    return vec

def spiral(x,y,z,k):
    A=1
    bas1=genvec([1,1,-2])
    bas2=genvec([1,1,1])
    r=N.array([x,y,z],'float64')
    u = A*(bas1[0]*N.cos(N.dot(k,r))+bas2[0]*N.sin(N.dot(k,r)))
    v = A*(bas1[1]*N.cos(N.dot(k,r))+bas2[1]*N.sin(N.dot(k,r)))
    w = A*(bas1[2]*N.cos(N.dot(k,r))+bas2[2]*N.sin(N.dot(k,r)))
    M=N.array([u,v,w])
    return M.T

def calc_struct(M,Qc,Q,r,P):
    """h,k,l are rhombohedral, hh,kh,lh are hexagonal"""
    s1=M

    flist=[]
    #find fperpendicular from F-F.q*q hat
    for i in range(len(Q)):
        #print 'i',i
        f=0
        vec=[Qc[i][0],Qc[i][1],Qc[i][2]]
        dotp1=N.dot(r,vec)*2*pi
        sp=s1-N.dot(s1,vec)*N.array(vec)/q[i]
        f=sp1*N.exp(-1.0j*dotp1)       
        #f2=N.abs(N.dot(f,N.conjugate(f)))

        
        #f2=f2/3
        #print 'f2inal',f2
        flist.append(f2)
    flist=N.array(flist)
    
    ff=mgnfacFe3psquared(q/4/pi)
    #print 'ff',ff
    Bf=.005
    flist=flist*ff*N.exp(-Bf*q**2/4/pi**2)
        
        
    return flist



def gen_atoms(r0,k,n):
    r=[]
    for i in range(-n,n):
        r.append(r0+k*i)
    r=N.array(r)
    return r
    
if __name__=="__main__":
    Qc=N.array([[0.5,0.5,0.5]])
    Q=setup(Qc[:,0],Qc[:,1],Qc[:,2])
    k=genvec([1,-1,0])*.0045
    r0=N.array([.5, .5, .5])
    n=100
    r=gen_atoms(r0,k,n)
    x=r[:,0]
    y=r[:,1]
    z=r[:,2]
    M=spiral(x,y,z,k)
    Pv=genvec(N.array([1,-1,0],'Float64'))
    Ph=genvec(N.array([1,1, 1],'Float64'))
    res=calc_struct(M,Qc,Q,r,Ph)
    print(M.shape)
    lam=2.35916


    
    
    