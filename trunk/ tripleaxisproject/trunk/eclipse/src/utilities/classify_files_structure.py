import readncnr3 as readncnr
import numpy as N
import scriptutil as SU
import re
import simple_combine



def readfiles(flist):
    mydatareader=readncnr.datareader()
    H=[]#N.array([])
    I=[]#N.array([])
    Ierr=[]#N.array([])
    monlist=[]
    count=0
    myfirstdata=mydatareader.readbuffer(flist[0])
    mon0=mydata.metadata['count_info']['monitor']
    print 'mon0',mon0
    
    for currfile in flist:
        print currfile
        mydata=mydatareader.readbuffer(currfile)


if __name__=='__main__':
    myfilestr=r'C:\Ce2RhIn8\Mar10_2009\magsc035.bt9'
    #myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\fpx53418.bt7'
    #myfilestr=r'c:\13165\13165\data\MagHigh56784.bt7'
    #myfilestr=r'c:\13176\data\CeOFeAs57255.bt7.out'
    mydatareader=readncnr.datareader()
    mydata=mydatareader.readbuffer(myfilestr)
    print mydata.__dict__
    print mydata.additional_metadata
    print mydata.metadata
    print mydata.metadata['file_info']['scantype']
    print mydata.metadata['collimations']
    print mydata.metadata['dspacing']['monochromator_dspacing']
    print mydata.metadata['dspacing']['analyzer_dspacing']
    print mydata.metadata['lattice']['a']
    print mydata.metadata['lattice']['b']
    print mydata.metadata['lattice']['c']
    print mydata.metadata['lattice']['alpha']
    print mydata.metadata['lattice']['beta']
    print mydata.metadata['lattice']['gamma']
    print mydata.metadata['motor3']['step']
    print mydata.metadata['motor4']['step']
    print mydata.metadata['count_info']['monitor']
    print mydata.metadata['energy_info']
    print mydata.metadata['q_center']['h_center']
    print mydata.metadata['q_center']['k_center']
    print mydata.metadata['q_center']['l_center']
    myfilebase='magsc'
    myend='bt9'
    mydirectory=r'c:\ce2rhin8\mar10_2009'
    myfilebaseglob=myfilebase+'*.'+myend
    print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    SU.printr(flist)
    
    