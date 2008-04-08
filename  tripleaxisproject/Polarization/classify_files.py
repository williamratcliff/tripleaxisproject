import numpy as N
import scriptutil as SU
import re
import readncnr2 as readncnr
import sys,os
threshold=1.0e-1

def readfiles(mydirectory,myfilebase,myend):
    myfilebaseglob=myfilebase+'*.'+myend
    #print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    mydatareader=readncnr.datareader()
    catalog={}
    catalog['pp']=[]
    catalog['mm']=[]
    catalog['pm']=[]
    catalog['mp']=[]
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
        data['filename']=currfile
        data['fully_polarized']=fully_polarized
        data['hsample']=hsample
        data['vsample']=vsample
        catalog[key].append(data)
        #print currfile, key
    return catalog



if __name__=='__main__':
    myend='bt7'
    mydirectory=r'c:\bifeo3xtal\jan8_2008\9175'
    myfilebase=''
    catalog=readfiles(mydirectory,myfilebase,myend)