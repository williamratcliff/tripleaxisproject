import numpy as N
import pylab
import scriptutil as SU
import re
import readncnr2 as readncnr
import simple_combine
#import scipy
from scipy.optimize import leastsq  
import copy
import scipy.odr
pi=N.pi

psd_a4_spacing=N.array([46.35,	  46.45,    46.54,     46.69,    46.84,     46.94,\
                        47.10,    47.23,    47.33,     47.46,    47.60,     47.71, \
                        47.85,    47.97,    48.14,     48.26,    48.33,     48.47, \
                        48.63,    48.72,    48.84,     48.95,    49.11,     49.21, \
                        49.38,    49.49,    49.61,     49.75,    49.90,     50.02, \
                        50.18,    50.28,    50.40,     50.53,    50.65,     50.81, \
                        50.91,    51.06,    51.19,     51.30,    51.46,     51.60, \
                        51.71,    51.83,    51.95,     52.09,    52.23,     52.14])
psd_a4_offset=psd_a4_spacing-psd_a4_spacing[22]
psd_a4_offset=N.reshape(psd_a4_offset,(psd_a4_offset.shape[0],1))
psd_channel_efficiency=N.array([3.9918,	2.4668,	2.7027,	2.5016,	2.2414,	2.2440,	2.0113,	1.900,\
                                1.7849,	1.7825,	1.6490,	1.4734,	1.3835,	1.3968,	1.3708,\
                                1.1635,	1.1048,	1.1963,	1.1978,	1.0392,	1.0331,	1.0618,\
                                1.0604,	1.0070,	0.9740,	1.0298,	1.0495,	0.9566,	0.9487,\
                                1.0109,	1.0164,	0.9770,	0.9357,	1.0383,	1.0691,	0.9457,\
                                0.9757,	0.9992,	1.0464,	0.9587,	0.9375,	0.9888,	0.9721,\
                                0.9391,	0.9335,	0.9982,	0.8798,	0.9205])

def read_order_files(flist):
#    myfilebaseglob=myfilebase+'*.'+myend
#    print myfilebaseglob
#    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    mydatareader=readncnr.datareader()
    H=[]#N.array([])
    I=[]#N.array([])
    Ierr=[]#N.array([])
    monlist=[]
    count=0
    mon0=5.0e4
    for currfile in flist:
        #print currfile
        mydata=mydatareader.readbuffer(currfile)
        #if count==0:
        #    #mon0=mydata.header['count_info']['monitor']
        #    mon0=5.0e4
        mon=mydata.metadata['count_info']['monitor']
        #print count, mon0,mon
        #temp=N.concatenate((temp,N.array(mydata.data['temp'])))
        H.append(N.array(mydata.data['magfield'],'float64'))
        psd_detectors=mydata.metadata['count_info']['AnalyzerDetectorDevicesOfInterest'.lower()]
        psd_arr=[]
        a4_center=[]
        a4=N.array(mydata.data['a4'],'float64')
        for psd in psd_detectors:
            psd_arr.append(N.array(mydata.data[psd],'float64'))
            a4_center.append(a4)
        print 'loop done'
        psd_arr=N.array(psd_arr,'float64')
        a4_center_arr=N.array(a4_center)
        print psd_arr.shape
        print a4_center_arr.shape
        #print psd_a4_offset
        print N.tile(psd_a4_offset,(1,a4.shape[0])).shape
        a4_corrected=a4_center_arr+N.tile(psd_a4_offset,(1,a4.shape[0]))
        print a4_corrected
        It=N.array(N.array(mydata.data['detector']),'float64')
        Iterr=N.sqrt(It)
        #It=It*mon0/mon
        monlist.append(mon)
        I.append(It)
        Ierr.append(Iterr)
        #Iterr=Iterr*mon0/mon
        #I=N.concatenate((I,It))
        #Ierr=N.concatenate((Ierr,Iterr))
        #print I
        #print Iterr
    #xa,ya,za=prep_data2(Qx,Qy,Counts);
        count=count+1
    return H,I,Ierr,monlist,a4


def orderparameter(p,T):
    #p[0]=Intensity
    I0,Tc,Beta,background=p
    I=I0*N.power(N.absolute(T/Tc-1),2*Beta)
    I[T>Tc]=0.0
    I=I+background
    return I

def gaussian(p,x):
    #p[0]=Intensity
    area,xc,w,background=p
    amp1=area/N.sqrt(2*pi)/w
    #[I Center Width Bak]
    peakshape=background+N.absolute(amp1)*N.exp(-0.5*((x-xc)/w)**2)
    return peakshape


def chisq_calc(p,x,I,Ierr):
    Icalc=gaussian(p,x)
    chisq=(I-Icalc)*(I-Icalc)/p.shape[0]/Ierr/Ierr
    return chisq

def residuals(p,x,I,Ierr):
    Icalc=gaussian(p,T)
    residual=(I-Icalc)/Ierr
    return residual

if __name__=='__main__':
    if 0:
        x=N.arange(-2,2,.25)
        p=[1.,0.5,1.,10.]
        y=gaussian(p,x)
        pylab.plot(x,y,linewidth=0,marker='s')
        pylab.show()
    if 1:
        mydirectory=r'C:\ca3comno6\Feb4_2008\Ca3CoMnO6\Feb4_2008\data'
        myfilebase='peak'
        myfileseq=54504
        myend='bt7'
        myfilestr=mydirectory+'\\'+myfilebase+str(myfileseq)+'.'+myend
        flist=[myfilestr]
        H,I,Ierr,monlist=read_order_files(flist)
        
        
    if 0:
        mydirectory=r'c:\camn2sb2\bt9\Feb5_2008'
        myfilebase='order*'
        myend='bt7'
        temp,I,Ierr,monlist=read_order_files(mydirectory,myfilebase,myend)
        T,I,Ierr=simple_combine.simple_combine(temp,I,Ierr,monlist)
        #print temp.shape, I.shape, Ierr.shape
        #pylab.errorbar(temp,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        p0=[copy.deepcopy(I[0]),84.0,.33333,copy.deepcopy(I[-1])]
        tmin=30
        tmax=100
        Trange=N.intersect1d(N.where(T>tmin)[0],N.where(T<tmax)[0])
        oparam=scipy.odr.Model(orderparameter)
        mydata=scipy.odr.RealData(T[Trange],I[Trange],sx=None,sy=Ierr)
        myodr = scipy.odr.ODR(mydata, oparam, beta0=p0)
        myoutput=myodr.run()
        myoutput.pprint()
        pfit=myoutput.beta
        #pfit = leastsq(residuals, p0, args=(T[Trange],I[Trange],Ierr[Trange]))
        print 'pfit=',pfit
        print 'chisq=',chisq_calc(pfit,T[Trange],I[Trange],Ierr[Trange]).sum()
        Icalc=orderparameter(pfit,T)
        if 1:
            pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
            pylab.plot(T,Icalc)
            pylab.xlabel('T (K)')
            pylab.ylabel('Counts (arb. units)')
            pylab.ylim((3000*0,7000))
            pylab.xlim((10,tmax))
            pylab.arrow(tmax,2000,0,500,fc='black',ec='black',width=.5)
            pylab.arrow(tmin,2000,0,500,fc='black',ec='black',width=.5)
            pylab.show()
