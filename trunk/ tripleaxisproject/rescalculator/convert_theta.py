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
import sys
sys.path.append(r'c:\tripleaxisproject2\polarization')
import polcorrect_lib
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
        data={}
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
        counts={}
        errs={}
        counts['pm']=I[0]
        counts['mp']=I[1]
        errs['pm']=Ierr[0]
        errs['mp']=Ierr[1]
        cell=mydirectory+'\WilliamOct2007horizCells.txt'
        print cell
        pbflags=polcorrect_lib.PBflags()
        pbflags.MonitorCorrect=0
        pbflags.PolMonitorCorrect=1
        pbflags.MonoSelect=1
        pbflags.Debug=0
        pbflags.SimFlux=0
        pbflags.SimDeviate=0
        pbflags.NoNegativeCS=0
        pbflags.HalfPolarized=0
        pbflags.CountsEnable=[0,0,1,1]
        pbflags.CountsAdd1=[0,0,0,0]
        pbflags.CountsAdd2=[0,0,0,0]
        pbflags.Sconstrain=[1,1,0,0]
        pbflags.Spp=[0,0,0,0]
        pbflags.Smm=[0,0,0,0]
        pbflags.Spm=[0,0,0,0]
        pbflags.Smp=[0,0,0,0]

        mypolcor=polcorrect_lib.polarization_correct(counts,errs,timestamp,cell,Ei_mp,Ef_mp)
        corrected_counts=mypolcor.correct(pbflags)

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



# Next zone
        myend='out'
        mydirectory=r'c:\bifeo3xtal\jan8_2008\9175\data'
        myfilebase='fieldscansplusminusreset53630'
        myfilebase2='fieldscanminusplusreset53631'

        myfilebaseglob=myfilebase+'*.'+myend
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilebaseglob=myfilebase2+'*.'+myend
        flist2 = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))

        myfilestr=flist[0]
        mydatareader=readncnr.datareader()
        mydata_pm=mydatareader.readbuffer(myfilestr)
        q_pm=N.array(mydata_pm.data['qx'])
        counts_pm=N.array(mydata_pm.data['detector_corrected'])
        errs_pm=N.array(mydata_pm.data['detector_errs_corrected'])


        myfilestr=flist2[0]
        mydatareader=readncnr.datareader()
        mydata_mp=mydatareader.readbuffer(myfilestr)
        q_mp=N.array(mydata_mp.data['qx'])
        counts_mp=N.array(mydata_mp.data['detector_corrected'])
        errs_mp=N.array(mydata_mp.data['detector_errs_corrected'])

        if 1:
            pylab.errorbar(q_mp,counts_mp,fmt='bo',yerr=errs_mp,linestyle='None')
            pylab.errorbar(q_pm,counts_pm,fmt='ro',yerr=errs_pm,linestyle='None')
            pylab.show()


        if 0:
            pylab.errorbar(H_mp,corrected_counts['Smp'],fmt='bo',yerr=corrected_counts['Emp'],linestyle='None')
            pylab.errorbar(H_pm,corrected_counts['Spm'],fmt='ro',yerr=corrected_counts['Epm'],linestyle='None')
            pylab.show()

        if 0:
            pylab.errorbar(H_mp,data['mp'],fmt='bo',yerr=errs['mp'],linestyle='None')
            pylab.errorbar(H_pm,data['pm'],fmt='ro',yerr=errs['pm'],linestyle='None')
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

