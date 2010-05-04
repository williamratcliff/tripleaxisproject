from utilities import readncnr3 as readncnr
import sys, os
import pylab
import numpy as N




if __name__=='__main__':
    mydir=r'c:\cepd\9072\data'
    myfilebase='CePd356111.bt7.out'
    myfilebase='CePd356123.bt7.out'
    myfilebase='CePd356141.bt7.out'
    myfilebase='CePd356142.bt7.out'
    myfilebase='CePd356144.bt7.out'
    myfilebase='CePd356145.bt7.out'
    myfilebase='CePd356111.bt7.out'
    myfilebase='CePd356111.bt7.out'
    myfilebase='CePd356129.bt7.out'
    myfilebase='CePd356138.bt7.out'
    #myfilebase='CePd356123.bt7.out'
    #myfilebase='CePd356149.bt7.out'
    myfilebase='CePd356151.bt7.out'
    myfilestr1=os.path.join(mydir,myfilebase)
    myfilebase='CePd356136.bt7.out'
    myfilebase='CePd356130.bt7.out'
    myfilebase='CePd356152.bt7.out'
    myfilebase='CePd356153.bt7.out'
    myfilebase='CePd356155.bt7.out'
    myfilebase='CePd356156.bt7.out'
    myfilebase='CePd356129.bt7.out'
    myfilebase='CePd356136.bt7.out'
    myfilebase='CePd356136.bt7.out'
    myfilebase='CePd356141.bt7.out'
    #myfilebase='CePd356130.bt7.out'
    #myfilebase='CePd356155.bt7.out'
    myfilebase='CePd356157.bt7.out'
    myfilestr2=os.path.join(mydir,myfilebase)   
    myreader=readncnr.datareader()
    mydata=myreader.readbuffer(myfilestr1)
    myreader=readncnr.datareader()
    mydata2=myreader.readbuffer(myfilestr2)
    #print 'hi'
    print mydata.metadata['count_info'].keys()
    signal1=mydata.metadata['count_info']['signal']
    signal2=mydata2.metadata['count_info']['signal']
    print signal1    
    varying1=mydata2.metadata['count_info']['varying'][0]
    varying2=mydata2.metadata['count_info']['varying'][0]
    print varying1

    mon1=N.array(mydata.data['monitor'])    
    mon2=N.array(mydata2.data['monitor'])
    
    mon0=mon1[0]
    factor1=mon0/mon1[0]
    factor2=mon0/mon2[0]
    print 'factor',factor2
    e1=N.array(mydata.data[varying1])
    e2=N.array(mydata2.data[varying2])
    I1=N.array(mydata.data[signal1+'_corrected'])*factor1
    I1err=N.array(mydata.data[signal1+'_errs_corrected'])*factor1
    
    I2=N.array(mydata2.data[signal2+'_corrected'])*factor2#*4/5
    I2err=N.array(mydata2.data[signal2+'_errs_corrected'])*factor2#*4/5
    
    #print 'E1',e1
    #print 'I',I1
    #print 'I1err',I1err
    #print 'E2',e2
    #print 'I2',I2
    #print 'I2err',I2err 
    print 'mon1', mon1 
    print 'mon2', mon2
    if 1:
        pylab.errorbar(e1,I1,I1err,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        pylab.errorbar(e2,I2,I2err,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
#    pylab.plot(e2,I2,'s','red')
        pylab.show()
    