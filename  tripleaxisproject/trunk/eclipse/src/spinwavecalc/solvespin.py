import numpy as N
import scipy
import math
import scipy.optimize
import scipy.optimize.anneal as anneal



def chisq_quaternion(p,a,b,c,a2,b2,c2,x0,y0,z0,x0p,y0p,z0p):
#For the quaternion case, the equations to be solved
    nq=N.sqrt(p[0]**2+p[1]**2+p[2]**2+p[3]**2)
    xq=p[0]
    yq=p[1]
    zq=p[2]
    wq=p[3]   
    cx=p[4]
    cy=p[5]
    cz=p[6]
    rmatrix=N.matrix([[1-2*yq**2-2*zq**2,2*xq*yq-2*wq*zq,2*zq*xq+2*wq*yq],\
                      [2*xq*yq+2*wq*zq,1-2*xq**2-2*zq**2,2*yq*zq-2*wq*xq],\
                      [2*zq*xq-2*wq*yq,2*yq*zq+2*wq*xq,1-2*xq**2-2*yq**2]],'d')
    translation=N.matrix([cx,cy,cz],'d').T
    xvec0=N.matrix([x0,y0,z0],'d').T
    xvec0p=transform_point(rmatrix,translation,xvec0)
    eqn4=normal_xq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn5=normal_yq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn6=normal_zq(xq,yq,zq,wq,cx,cy,cz,a,b,c,a2,b2,c2,x0,y0,z0)
    eqn1=N.ravel(xvec0p[0]-x0p)[0]
    eqn2=N.ravel(xvec0p[1]-y0p)[0]
    eqn3=N.ravel(xvec0p[2]-z0p)[0]
    normq=xq**2+yq**2+zq**2+wq**2-1  #require a unit quaternion
    fresult=N.array([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,normq],'d')
    return fresult

if __name__=="__main__":
    p0=N.array([0,-1,0,0,6,0,0],'d')

#use a simple nonlinear solver for finding A solution to our system of equations 
    p=scipy.optimize.minpack.fsolve(chisq_quaternion,p0,args=(a,b,c,a2,b2,c2,x0,y0,z0,x0p,y0p,z0p))
