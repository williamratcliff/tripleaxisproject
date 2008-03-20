import numpy as N
import pylab
import math
import unittest
from matplotlib.patches import Ellipse
import lattice_calculator
import readncnr2 as readncnr
import scriptutil as SU
import re
eps=1e-3
pi=N.pi


def autovectorized(f):
     """Function decorator to do vectorization only as necessary.
     vectorized functions fail for scalar inputs."""
     def wrapper(input):
         if N.isscalar(input)==False:
             return N.vectorize(f)(input)
         return f(input)
     return wrapper



@autovectorized
def myradians(x):
    return math.radians(x)

def readfiles(mydirectory,myfilebase,myend):
    myfilebaseglob=myfilebase+'*.'+myend
    print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))


if __name__=="__main__":
    if 1:
        a=N.array([5.5436],'d')
        b=N.array([5.5436],'d')
        c=N.array([13.8844],'d')
        alpha=N.array([90],'d')
        beta=N.array([90],'d')
        gamma=N.array([120],'d')
 #       orient1=N.array([[0,1,1]],'d')
        orient1=N.array([[0,0,1]],'d')
        orient2=N.array([[1,1,0]],'d')
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)

        mydirectory=r'c:\BiFeO3xtal\Oct9_2007'
        myfilebase='bfo_spinflip_he3flip51585'
        myfilebase2='bfospinflip51579'
        myfilebase3='bfospinflip51583'
        myend='bt7'
        myfilebaseglob=myfilebase+'*.'+myend
        print myfilebaseglob
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilestr=flist[0]
        Counts=N.array([],'float64')
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        S1=N.array(mydata.data['a3'])
        S2=N.array(mydata.data['a4'])
        A2=N.array(mydata.data['a5'])*2
        M2=N.array(mydata.data['a2'])
        Counts=N.concatenate((Counts,N.array(mydata.data['detector'])))
        if 0:
            for i in range(S1.size):
                s='%f %f %f %f'%(S1[i],S2[i],A2[i],M2[i])
                print s
            exit()
        if 0:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians([M2[0]]),myradians([-163.193]),myradians([29.833]),myradians([A2[0]]),setup)
            s='%f %f %f %f %f %f %f'%(S1[0],S2[0],H[0],K[0],L[0],E[0],Q[0])
            print s
        if 1:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
            for i in range(1,S1.size):
                s='%f %f %f %f %f %f %f'%(S1[i],S2[i],H[i],K[i],L[i],E[i],Q[i])
                print s

        if 1:
            pylab.plot(H,L,'bo')
            pylab.xlabel('[1,1,0]')
            pylab.ylabel('[0,0,1]')
            #pylab.axis([-1,1,-4,4])
            pylab.show()

