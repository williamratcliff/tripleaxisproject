import numpy as N
import scriptutil as SU
import re
import readncnr2 as readncnr
import sys,os
threshold=1.0e-1

#plus refers to a flipper state that is on
#minus refers to a flipper state that is off

class pol_info:
    def __init__(self):
        self.data=[]
        self.files=[]
        return


class catalog:
    def __init__(self):
        self.pp=pol_info()
        self.pm=pol_info()
        self.mp=pol_info()
        self.mm=pol_info()
        return

def readfiles(flist):
    #SU.printr(flist)
    mydatareader=readncnr.datareader()
    mycatalog=catalog()
    half_polarized=0
    for currfile in flist:
        #print currfile
        key_i='m'
        hsample='none'
        vsample='none'
        mydata=mydatareader.readbuffer(currfile,lines=3)
        if mydata.data.has_key('hsample'):
            hsample=mydata.data['hsample'][0]
        if mydata.data.has_key('vsample'):
            vsample=mydata.data['vsample'][0]
        ei_exists=0
        if mydata.data.has_key('eiflip'):
            ei_exists=1
            #print mydata.data['eiflip']
            if N.average(mydata.data['eiflip'])>threshold:
                key_i='p'
        key_f='m'
        ef_exists=0
        if mydata.data.has_key('efflip'):
            ef_exists=1
            if N.average(mydata.data['efflip'])>threshold:
                key_f='p'
        key=key_i+key_f
        fully_polarized=ei_exists*ef_exists  #is 1 if both cells are present
        data={}
        data['absolute_filename']=currfile
        data['fully_polarized']=fully_polarized
        data['hsample']=hsample
        data['vsample']=vsample
        pattern = re.compile('^(?P<base>[^.]*?)(?P<seq>[0-9]*)(?P<ext>[.].*)?$')
        filename=os.path.split(currfile)[-1]
        #print filename
        match = pattern.match(filename)
        dict((a,match.group(a)+"") for a in ['base','seq','ext'])
        data['filebase']=match.group('base')
        data['fileseq_number']=match.group('seq')
        data['filename']=data['filebase']+data['fileseq_number']

        if key=='pp':
            mycatalog.pp.data.append(data)
            mycatalog.pp.files.append(data['filename'])
        elif key=='pm':
            mycatalog.pm.data.append(data)
            mycatalog.pm.files.append(data['filename'])
        elif key=='mp':
            mycatalog.mp.data.append(data)
            mycatalog.mp.files.append(data['filename'])
        elif key=='mm':
            mycatalog.mm.data.append(data)
            mycatalog.mm.files.append(data['filename'])

        #print currfile, key
    return mycatalog



if __name__=='__main__':
    myend='bt7'
    mydirectory=r'c:\bifeo3xtal\jan8_2008\9175'
    myfilebase=''
    myfilebaseglob=myfilebase+'*.'+myend
    #print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    mycatalog=readfiles(flist)
    print mycatalog.pm.files