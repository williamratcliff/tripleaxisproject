import csv
import numpy as N
import array
import pylab

pi=N.pi

def fwhm2sigma(fwhm):
    factor=2*N.sqrt(2*N.log(2.0))
    sigma=fwhm/factor
    return sigma

def area2I(A,fwhm):
    sigma=fwhm2sigma(fwhm)
    print 'sigma ',sigma
    I=A/(N.sqrt(2*pi)*sigma)
    return I

def readfile(myfilestr):
    reader = csv.reader(open(myfilestr))
    header = reader.next()
    data = array.array('f')
    for row in reader:
        data.extend(map(float, row))
    #print 'Data size', len(data)
    #myfile=open(myfilestr)
    #header = myfile.readline()
    #data = N.fromfile(myfile, sep=',')#.reshape(-1,3)
    #myfile.close()
    ndata=N.array(data).reshape(-1,3)
    return ndata

def calcrange(a4lim,data):
    """returns the range of the data array included in the a4lim tuple (a4min,a4max)"""
    a4range=N.intersect1d(N.where(data>a4lim[0])[0],N.where(data<a4lim[1])[0])
    return a4range



if __name__=='__main__':

    print area2I(87.82,.4537)
    mydirectory=r'c:\\12436\superconductor'
    myend='csv'
    myfilebase='lafen015_220_'
    T=[4,134,138,141,148,153,155,175]
    a4lim=(92.0,95.0)
    a4=N.array([],'float64')
    I=N.array([],'float64')
    Ierr=N.array([],'float64')
    Tlist=N.array([],'float64')
    for currT in T:
        myfilestr=mydirectory+'\\'+str(currT)+'K.'+myend
        print myfilestr
        data=readfile(myfilestr)
        if currT==4:
            data[:,1]=data[:,1]/3
            data[:,2]=data[:,2]/3
            #print '4'
        if currT==60:
            data[:,1]=data[:,1]/6
            data[:,2]=data[:,2]/6
            print 60
        if currT==175:
            data[:,1]=data[:,1]*0.25
            data[:,2]=data[:,2]*0.25
            #print '175'
        #data[:,0]->a4, others go [a4,I,Ierr]
        a4range=calcrange(a4lim,data[:,0])
        if 1:
            if currT in[175,4]:
                pylab.errorbar(data[a4range,0],data[a4range,1],data[a4range,2],marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        a4=N.concatenate((a4,data[a4range,0]))
        I=N.concatenate((I,data[a4range,1]))
        Ierr=N.concatenate((Ierr,data[a4range,2]))
        Tlist=N.concatenate((Tlist,currT*N.ones(data[a4range,2].shape,'float64')))
    if 0:
        pylab.show()
    print a4
    a4=N.array(a4).flatten()
    I=N.array(I).flatten()
    Ierr=N.array(Ierr).flatten()
    T=N.array(Tlist).flatten()

    print T

    if 0:
        pylab.errorbar(data[a4range,0],data[a4range,1],data[a4range,2],marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        pylab.show()