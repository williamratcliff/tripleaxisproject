import numpy as N
import pylab
import math
import unittest
from matplotlib.patches import Ellipse
import lattice_calculator
import readncnr2 as readncnr
import scriptutil as SU
import re
import simple_combine
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
    Qx=N.array([])
    Qy=N.array([])
    Qz=N.array([])
    Counts=N.array([])
    for currfile in flist:
        print currfile
        mydata=mydatareader.readbuffer(currfile)
        Qx=N.concatenate((Qx,N.array(mydata.data['qx'])))
        Qy=N.concatenate((Qy,N.array(mydata.data['qy'])))
        Qz=N.concatenate((Qz,N.array(mydata.data['qz'])))
        Counts=N.concatenate((Counts,N.array(mydata.data['counts'])))


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

        EXP={}
        EXP['ana']={}
        EXP['ana']['tau']='pg(002)'
        EXP['mono']={}
        EXP['mono']['tau']='pg(002)';
        EXP['ana']['mosaic']=25
        EXP['mono']['mosaic']=25
        EXP['sample']={}
        EXP['sample']['mosaic']=25
        EXP['sample']['vmosaic']=25
        EXP['hcol']=N.array([40, 40, 40, 80],'d')
        EXP['vcol']=N.array([120, 120, 120, 120],'d')
        EXP['infix']=-1 #positive for fixed incident energy
        EXP['efixed']=14.7
        EXP['method']=0
        setup=[EXP]

        mydirectory=r'c:\BiFeO3xtal\Oct9_2007'
        myfilebase='bfo_spinflip_he3flip51585'
        myfilebase2='bfo_spinflip51579'
        myfilebase3='bfo_spinflip51583'
        myend='bt7'

        #pm
        myfilebaseglob=myfilebase+'*.'+myend
        print myfilebaseglob
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilestr=flist[0]
        Counts=N.array([],'float64')
        monlist=[]
        timestamp={}
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        S1=N.array(mydata.data['a3'][1:])
        S2=N.array(mydata.data['a4'][1:])
        A2=N.array(mydata.data['a5'][1:])*2
        M2=N.array(mydata.data['a2'][1:])
        Counts_pm=N.concatenate((Counts,N.array(mydata.data['detector'][1:])))
        H_pm,K_pm,L_pm,E_pm,Q_pm,Ei_pm,Ef_pm=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
        monlist.append(mydata.data['monitor'][0])
        Counts_pm=Counts_pm[:-1]
        H_pm=H_pm[:-1]
        K_pm=K_pm[:-1]
        L_pm=L_pm[:-1]
        E_pm=E_pm[:-1]
        Q_pm=Q_pm[:-1]
        Ei_pm=Ei_pm[:-1]
        Ef_pm=Ef_pm[:-1]
        timestamp['pm']=N.array(mydata.data['timestamp'][1:])[:-1]
        #print H_pm
        #begin mp
        myfilebaseglob=myfilebase2+'*.'+myend
        print myfilebaseglob
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilebaseglob=myfilebase3+'*.'+myend
        flist2 = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))

        myfilestr=flist[0]
        Counts=N.array([],'float64')
        mydatareader=readncnr.datareader()
        mydata1=mydatareader.readbuffer(myfilestr)
        S1=N.array(mydata1.data['a3'])[1:]
        S2=N.array(mydata1.data['a4'])[1:]
        A2=N.array(mydata1.data['a5'])[1:]*2
        M2=N.array(mydata1.data['a2'])[1:]
        Counts_mp=N.array(mydata1.data['detector'])[1:]
        timestamp_mp=N.array(mydata1.data['timestamp'])[1:]
        #H_mp,K_mp,L_mp,E_mp,Q_mp,Ei_mp,Ef_mp=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
        monlist.append(mydata1.data['monitor'][0])
        #print M2.shape
        #print Counts_mp.shape
        #remove potential bad datapoint that bt7 generated in dd mode


        myfilestr=flist2[0]
        mydatareader=readncnr.datareader()
        mydata2=mydatareader.readbuffer(myfilestr)
        S1=N.concatenate((S1,N.array(mydata2.data['a3'])[1:]))
        S2=N.concatenate((S2,N.array(mydata2.data['a4'])[1:]))
        A2=N.concatenate((A2,N.array(mydata2.data['a5'])[1:]*2))
        M2=N.concatenate((M2,N.array(mydata2.data['a2'])[1:]))
        timestamp_mp=N.concatenate((timestamp_mp,N.array(mydata2.data['timestamp'])[1:]))
        timestamp['mp']=timestamp_mp
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)

        Counts_mp=N.concatenate((Counts_mp,N.array(mydata2.data['detector'])[1:]))
        #print 'where'
        H_mp,K_mp,L_mp,E_mp,Q_mp,Ei_mp,Ef_mp=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)


        ylist=[Counts_pm,Counts_mp]
        yerrlist=[N.sqrt(Counts_pm),N.sqrt(Counts_mp)]
        I,Ierr=simple_combine.monitor_normalize(ylist,yerrlist,monlist)
        data={}
        data['pm']=I[0]
        data['mp']=I[1]
        #print timestamp['mp'].shape
        #print timestamp['pm'].shape
        #print S1.shape
        #print S2.shape
        #print A2.shape
        #print M2.shape
        #print Counts_mp.shape
        #print H_mp.shape
        #print Counts_mp.shape
        #print I[0].shape
        #print I[1].shape
        if 0:
            pylab.plot(H_mp,data['mp'],'bo')
            pylab.plot(H_pm,data['pm'],'ro')
            pylab.show()


        exit()
        if 0:
            for i in range(S1.size):
                s='%f %f %f %f'%(S1[i],S2[i],A2[i],M2[i])
                print s
            exit()
        if 0:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians([M2[0]]),myradians([-163.193]),myradians([29.833]),myradians([A2[0]]),setup)
            s='%f %f %f %f %f %f %f'%(S1[0],S2[0],H[0],K[0],L[0],E[0],Q[0])
            print s
        if 0:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
            for i in range(1,S1.size):
                s='%f %f %f %f %f %f %f'%(S1[i],S2[i],H[i],K[i],L[i],E[i],Q[i])
                print s

        if 0:
            pylab.plot(H,L,'bo')
            pylab.xlabel('[1,1,0]')
            pylab.ylabel('[0,0,1]')
            #pylab.axis([-1,1,-4,4])
            pylab.show()

