from datetime import datetime
import time
import sys,os
import numpy as np
import copy
import ctypes
import cstruct
import os, sys
from ctypes import c_void_p, c_int, c_long, c_char, c_char_p,c_double,c_ulong
from ctypes import byref as _ref
c_void_pp = ctypes.POINTER(c_void_p)
c_char_pp = ctypes.POINTER(c_char_p)
c_int_p = ctypes.POINTER(c_int)
c_long_p=ctypes.POINTER(c_long)
c_double_p=ctypes.POINTER(c_double)
c_ulong_p=ctypes.POINTER(c_ulong)



#CTYPES STUFF

##PBdatapt = cstruct.cstruct(('sec',N.float64,(4,1)),('Y',N.float64,(4,1)),
##                    ('Yesq',N.float64,(4,1)),('lambI',N.float64,(4,1)),
##                    ('lamF',N.float64,(4,1)),('C',N.float64,(4,4)),
##                    ('Cesq',N.float64,(4,4)),('S',N.float64,(4,1)),
##                    ('Sesq',N.float64,(4,1)),('NActive',N.float64,(4,1)),
##                    ('activeEq',N.float64,(4,1)),('Nfree',N.intc),
##                    ('freeS',numpy.intc))


class PBindata(ctypes.Structure):
    _fields_=[("Ei",c_double_p),
            ("Ef",c_double_p),
            ("Cpp",c_double_p),
            ("Cmm",c_double_p),
            ("Cpm",c_double_p),
            ("Cmp",c_double_p),
            ("Epp",c_double_p),
            ("Emm",c_double_p),
            ("Epm",c_double_p),
            ("Emp",c_double_p),
            ("tpp",c_ulong_p),
            ("tmm",c_ulong_p),
            ("tpm",c_ulong_p),
            ("tmp",c_ulong_p)
            ]

#typedef struct {
#  double *Ei, *Ef ;
#  double *Cpp, *Cmm, *Cpm, *Cmp ;
#  double *Epp, *Emm, *Epm, *Emp ;
#  unsigned long *tpp, *tmm, *tpm, *tmp ;
#} PBindata ;

class PBoutdata(ctypes.Structure):
    _fields_=[("Spp",c_double_p),
            ("Smm",c_double_p),
            ("Spm",c_double_p),
            ("Smp",c_double_p),
            ("Epp",c_double_p),
            ("Emm",c_double_p),
            ("Epm",c_double_p),
            ("Emp",c_double_p),
            ("R",c_double_p)
            ]

#typedef struct {
#  double *Spp, *Smm, *Spm, *Smp ;
#  double *Epp, *Emm, *Epm, *Emp ;
#} PBoutdata ;

##class PBflags(ctypes.Structure):
##    _fields_=[("MonitorCorrect",c_int),
##            ("PolMonitorCorrect",c_int),
##            ("Debug",c_int),
##            ("SimFlux",c_int),
##            ("SimDeviate",c_int),
##            ("CountsEnable",c_int_p),
##            ("CountsAdd1",c_int_p),
##            ("CountsAdd2",c_int_p),
##            ("Sconstrain",c_int_p),
##            ("Spp",c_double_p),
##            ("Smm",c_double_p),
##            ("Spm",c_double_p),
##            ("Smp",c_double_p)
##            ]

PBflags=cstruct.cstruct(('MonitorCorrect',N.intc),
                    ('PolMonitorCorrect',N.intc),
                    ('MonoSelect',N.intc),
                    ('Debug',N.intc),
                    ('SimFlux',N.intc),
                    ('SimDeviate',N.intc),
                    ('NoNegativeCS',N.intc),
                    ('HalfPolarized',N.intc),
                    ('CountsEnable',int,(4,)),
                    ('CountsAdd1',int,(4,)),
                    ('CountsAdd2',int,(4,)),
                    ('Sconstrain',int,(4,)),
                    ('Spp','float64',(4,)),
                    ('Smm','float64',(4,)),
                    ('Spm','float64',(4,)),
                    ('Smp','float64',(4,))
                    )




#typedef struct {
#  int MonitorCorrect ;
#  int PolMonitorCorrect ;
#  int Debug ;
#  int SimFlux ;
#  int SimDeviate ;
#  int CountsEnable[4] ;
#  int CountsAdd1[4] ;
#  int CountsAdd2[4] ;
#  int Sconstrain[4] ;
#  double Spp[4], Smm[4], Spm[4], Smp[4] ;
#} PBflags ;


#/* entrypoints */

#int PBcorrectData(char *PCellFile, char *ACellFile, PBflags *flgs,
#		  int npts, PBindata *in, PBoutdata *out) ;
#int PBsim(char *filename) ;
#int PBreadflags(char *filename) ;

if sys.platform=='win32':
    mypolcorrect = N.ctypeslib.load_library('polarization2.dll', 'C:\mytripleaxisproject\tripleaxisproject\eclipse\src\polarization\sanspol')
    #mypolcorrect = N.ctypeslib.load_library(r'C:\mytripleaxisproject\trunk\eclipse\src\polarization\polarization2.dll')
elif sys.platform=='mac':
    mypolcorrect = N.ctypeslib.load_library('libpolarization2.so', '.')
else:
    mypolcorrect = N.ctypeslib.load_library('libpolarization2.so', '.') #linux


    def correct(self,pbflags):
        #print pbflags.Sconstrain
        keys=['pp','mm','pm','mp']
        #print self.counts
        #print self.timestamp
        #pbinput=PBindata()
        #pboutput=PBoutdata()
        #pbinput.Ei=self.ei.ctypes.data_as(c_double_p)
        #binput.Ef=self.ef.ctypes.data_as(c_double_p)
        #print 'Ei', self.ei
        #print 'Ef',self.ef
                    

        self.outdata={}
        for key in keys:
            if self.counts.has_key(key):
                lastkey=key
                self.outdata[key]=copy.deepcopy(self.mydata[key])
                detectors=self.mydata[key].metadata['count_info']['AnalyzerDetectorDevicesOfInterest'.lower()]
        detectors.append(self.mydata[lastkey].metadata['count_info']['signal'])
        for detector in detectors: 
        #detector=detectors[-1]
        #if 1:
            pbinput=PBindata()
            pboutput=PBoutdata()
            pbinput.Ei=self.ei.ctypes.data_as(c_double_p)
            pbinput.Ef=self.ef.ctypes.data_as(c_double_p)
            pbflags.MonitorCorrect=int(pbflags.MonitorCorrect)
            pbflags.PolMonitorCorrect=int(pbflags.PolMonitorCorrect)
            pbflags.MonoSelect=int(pbflags.MonoSelect)
            pbflags.Debug=int(pbflags.Debug)
            pbflags.SimFlux=int(pbflags.SimFlux)
            pbflags.SimDeviate=int(pbflags.SimDeviate)
            pbflags.NoNegativeCS=int(pbflags.NoNegativeCS)
            pbflags.HalfPolarized=int(pbflags.HalfPolarized)
            pbflags.CountsEnable[0]=int(pbflags.CountsEnable[0])
            pbflags.CountsEnable[1]=int(pbflags.CountsEnable[1])
            pbflags.CountsEnable[2]=int(pbflags.CountsEnable[2])
            pbflags.CountsEnable[3]=int(pbflags.CountsEnable[3])
            pbflags.CountsAdd1[0]=int(pbflags.CountsAdd1[0])
            pbflags.CountsAdd1[1]=int(pbflags.CountsAdd1[1])
            pbflags.CountsAdd1[2]=int(pbflags.CountsAdd1[2])
            pbflags.CountsAdd1[3]=int(pbflags.CountsAdd1[3])
            pbflags.CountsAdd2[0]=int(pbflags.CountsAdd2[0])
            pbflags.CountsAdd2[1]=int(pbflags.CountsAdd2[1])
            pbflags.CountsAdd2[2]=int(pbflags.CountsAdd2[2])
            pbflags.CountsAdd2[3]=int(pbflags.CountsAdd2[3])
    
            pbflags.Sconstrain[0]=int(pbflags.Sconstrain[0])
            pbflags.Sconstrain[1]=int(pbflags.Sconstrain[1])
            pbflags.Sconstrain[2]=int(pbflags.Sconstrain[2])
            pbflags.Sconstrain[3]=int(pbflags.Sconstrain[3])
#            pbflags.Spp[0]=float(pbflags.Spp[0])
#            pbflags.Spp[1]=float(pbflags.Spp[1])
#            pbflags.Spp[2]=float(pbflags.Spp[2])
#            pbflags.Spp[3]=float(pbflags.Spp[3])
#            pbflags.Smm[0]=float(pbflags.Smm[0])
#            pbflags.Smm[1]=float(pbflags.Smm[1])
#            pbflags.Smm[2]=float(pbflags.Smm[2])
#            pbflags.Smm[3]=float(pbflags.Smm[3])
#            pbflags.Spm[0]=float(pbflags.Spm[0])
#            pbflags.Spm[1]=float(pbflags.Spm[1])
#            pbflags.Spm[2]=float(pbflags.Spm[2])
#            pbflags.Spm[3]=float(pbflags.Spm[3])
#            pbflags.Smp[0]=float(pbflags.Smp[0])
#            pbflags.Smp[1]=float(pbflags.Smp[1])
#            pbflags.Smp[2]=float(pbflags.Smp[2])
#            pbflags.Smp[3]=float(pbflags.Smp[3])
    
            if self.mydata[lastkey].metadata['count_info']['count_type']=='time':
                pbflags.MonitorCorrect=int(0)
                pbflags.PolMonitorCorrect=int(0)




 
            for key in keys:
                if self.counts.has_key(key):
                    mon0=self.mon0
                    mon=self.mon[key]
                    #print 'mon0',mon0
                    #print 'mon',mon                 
                    monfactor=N.array(mon0,'float64')/N.array(mon,'float64')
                    if self.mydata[lastkey].metadata['count_info']['count_type']=='time':
                        monfactor=1.0
                    #print 'monfactor', monfactor
                    self.counts[key]=N.array(self.mydata[key].data[detector],'float64')*monfactor
                    #print 'before norm'
                    #print N.array(self.mydata[key].data[detector])
                    self.errors[key]=N.sqrt(N.array(self.mydata[key].data[detector],'float64'))*monfactor 
                    #print 'counts',self.counts[key]
                    #print 'errors',self.errors[key] 
                    if key=='pp':
                        mytemp_pp=self.timestamp[key].astype('uint32')
                        pbinput.tpp=mytemp_pp.ctypes.data_as(c_ulong_p)
                        #pbinput.tpp=self.timestamp[key].astype('uint32').ctypes.data_as(c_ulong_p)
                        pbinput.Cpp=self.counts[key].ctypes.data_as(c_double_p)
                        pbinput.Epp=self.errors[key].ctypes.data_as(c_double_p)
                    if key=='mm':
                        mytemp_mm=self.timestamp[key].astype('uint32')
                        pbinput.tmm=mytemp_mm.ctypes.data_as(c_ulong_p)
                        #pbinput.tmm=self.timestamp[key].astype('uint32').ctypes.data_as(c_ulong_p)
                        #print 't0,mm ',mytemp_mm[0]
                        #print 'counts',self.counts[key][0]
                        pbinput.Cmm=self.counts[key].ctypes.data_as(c_double_p)
                        pbinput.Emm=self.errors[key].ctypes.data_as(c_double_p)
                    if key=='pm':
                        mytemp_pm=self.timestamp[key].astype('uint32')
                        pbinput.tpm=mytemp_pm.ctypes.data_as(c_ulong_p)
                        #print 't0,pm ',mytemp_pm[0]
                        #print 'counts',self.counts[key][0]
                        pbinput.Cpm=self.counts[key].ctypes.data_as(c_double_p)
                        pbinput.Epm=self.errors[key].ctypes.data_as(c_double_p)
                        #print 'shape ',self.counts[key].shape
                    if key=='mp':
                        mytemp_mp=self.timestamp[key].astype('uint32')
                        pbinput.tmp=mytemp_mp.ctypes.data_as(c_ulong_p)
                        pbinput.Cmp=self.counts[key].ctypes.data_as(c_double_p)
                        pbinput.Emp=self.errors[key].ctypes.data_as(c_double_p)    
                else:    
                    dummytpp=N.empty((1,self.length),'uint32')
                    dummytmm=N.empty((1,self.length),'uint32')
                    dummytpm=N.empty((1,self.length),'uint32')
                    dummytmp=N.empty((1,self.length),'uint32')
                    dummyCpp=N.empty((1,self.length),'float64')
                    dummyCmm=N.empty((1,self.length),'float64')
                    dummyCpm=N.empty((1,self.length),'float64')
                    dummyCmp=N.empty((1,self.length),'float64')
                    dummyEpp=N.empty((1,self.length),'float64')
                    dummyEmm=N.empty((1,self.length),'float64')
                    dummyEpm=N.empty((1,self.length),'float64')
                    dummyEmp=N.empty((1,self.length),'float64')
    
                    if key=='pp':
                        pbinput.tpp=dummytpp.ctypes.data_as(c_ulong_p)
                        pbinput.Cpp=dummyCpp.ctypes.data_as(c_double_p)
                        pbinput.Epp=dummyEpp.ctypes.data_as(c_double_p)
                    if key=='mm':
                        pbinput.tmm=dummytmm.ctypes.data_as(c_ulong_p)
                        pbinput.Cmm=dummyCmm.ctypes.data_as(c_double_p)
                        pbinput.Emm=dummyEmm.ctypes.data_as(c_double_p)
                    if key=='pm':
                        pbinput.tpm=dummytpm.ctypes.data_as(c_ulong_p)
                        pbinput.Cpm=dummyCpm.ctypes.data_as(c_double_p)
                        pbinput.Epm=dummyEpm.ctypes.data_as(c_double_p)
                    if key=='mp':
                        pbinput.tmp=dummytmp.ctypes.data_as(c_ulong_p)
                        pbinput.Cmp=dummyCmp.ctypes.data_as(c_double_p)
                        pbinput.Emp=dummyEmp.ctypes.data_as(c_double_p)
    
        
    
    
            Spp=N.empty((1,self.length),'float64')
            Smm=N.empty((1,self.length),'float64')
            Spm=N.empty((1,self.length),'float64')
            Smp=N.empty((1,self.length),'float64')
    
            pboutput.Spp=Spp.ctypes.data_as(c_double_p)
            pboutput.Smm=Smm.ctypes.data_as(c_double_p)
            pboutput.Spm=Spm.ctypes.data_as(c_double_p)
            pboutput.Smp=Smp.ctypes.data_as(c_double_p)
    
            Epp=N.empty((1,self.length),'float64')
            Emm=N.empty((1,self.length),'float64')
            Epm=N.empty((1,self.length),'float64')
            Emp=N.empty((1,self.length),'float64')
            R=N.empty((1,self.length),'float64')
    
            pboutput.Epp=Epp.ctypes.data_as(c_double_p)
            pboutput.Emm=Emm.ctypes.data_as(c_double_p)
            pboutput.Epm=Epm.ctypes.data_as(c_double_p)
            pboutput.Emp=Emp.ctypes.data_as(c_double_p)
            pboutput.R=R.ctypes.data_as(c_double_p)

            
#            pbflags=PBflags()
#    
#            fMonitorCorrect=0
#            fPolMonitorCorrect=0
#            if self.mydata[lastkey].metadata['count_info']['count_type']=='monitor':
#                fPolMonitorCorrect=1
#    ##            #print 'monitor'
#    ##        fDebug=0
#    ##        fSimFlux=0
#    ##        fSimDeviate=0
#            #-- ++ +- -+ #TODO Check with Ross on the order
#    
#    
#            pbflags.MonitorCorrect=0#fMonitorCorrect
#            pbflags.PolMonitorCorrect=0#fPolMonitorCorrect
#            pbflags.MonoSelect=1
#            pbflags.Debug=9#fDebug
#            pbflags.SimFlux=0#fSimFlux
#            pbflags.SimDeviate=0#fSimDeviate
#            pbflags.NoNegativeCS=0
#            pbflags.HalfPolarized=0
#            pbflags.CountsEnable=[0,1,1,0]
#            pbflags.CountsAdd1=[0,0,0,0]
#            pbflags.CountsAdd2=[0,0,0,0]
#            pbflags.Sconstrain=[1,0,0,1]
#            pbflags.Spp=[0,1,0,0]
#            pbflags.Smm=[0,0,0,0]
#            pbflags.Spm=[0,0,0,0]
#            pbflags.Smp=[0,0,1,0]
#            self.cell=r'c:\13176\data\CeOFeCellsMay20081.txt'
        
#            print 'Debug',pbflags.Debug
#            print 'SimFlux',pbflags.SimFlux
#            print 'SimDeviate',pbflags.SimDeviate
#            print 'NoNegativeCS',pbflags.NoNegativeCS
#            print 'HalfPolarized', pbflags.HalfPolarized
#            print 'CountsAdd1',pbflags.CountsAdd1
#            print 'CountsAdd2',pbflags.CountsAdd2
#            print 'monitor correct ', pbflags.MonitorCorrect
#            print 'mono select ',pbflags.MonoSelect
#            print 'PolMonitorCorrect',pbflags.PolMonitorCorrect
#            print 'Spm ',pbflags.Spm
#            print 'Smp ',pbflags.Smp
#            print 'Smm ',pbflags.Smm
#            print 'Spp ',pbflags.Spp
#            print 'Sconstrain ', pbflags.Sconstrain
#            print 'CountsEnable ', pbflags.CountsEnable
#            print 'mycell', self.cell
            ierr=mypolcorrect.PBcorrectData(self.cell,pbflags._pointer,self.length,ctypes.byref(pbinput),ctypes.byref(pboutput))
            #print Smm[0]
            #print Spm[0]
            #print 'ierr',ierr
            corrected_counts={}
            corrected_counts['Spp']=Spp[0]
            corrected_counts['Epp']=Epp[0]
    
            corrected_counts['Smm']=Smm[0]
            corrected_counts['Emm']=Emm[0]
    
            corrected_counts['Spm']=Spm[0]
            corrected_counts['Epm']=Epm[0]
    
            corrected_counts['Smp']=Smp[0]
            corrected_counts['Emp']=Emp[0]
            self.corrected_counts=corrected_counts

            for key in keys:
                if self.counts.has_key(key):
                    #self.outdata[key]=copy.deepcopy(self.mydata[key])
                    ##self.mydata[key].data[self.mydata[key].metadata['signal']]
                    newfield=detector+'_corrected'
                    #print 'key',key,'newfield',newfield
                    #print 'corrected counts',corrected_counts['S'+key]
                    #print 'corrected Errors',corrected_counts['E'+key]
                    self.outdata[key].data[newfield]=corrected_counts['S'+key]
                    newfield=detector+'_errs_corrected'
                    self.outdata[key].data[newfield]=corrected_counts['E'+key]
                    #print 'correcting key ',key,' len ',self.outdata[key].data[newfield].shape
            #corrected_counts_out=corrected_counts
            if detector==self.mydata[lastkey].metadata['count_info']['signal']:
                corrected_counts_out=corrected_counts
        return corrected_counts_out



def get_tokenized_line(myfile,returnline=[''],splitchar=None):
    lineStr=myfile.readline()
    returnline[0]=lineStr.rstrip()
    strippedLine=lineStr.lower().rstrip().lstrip()
    if splitchar==None:
        tokenized=strippedLine.split()
    else:
        tokenized=strippedLine.split(splitchar)
    return tokenized


def fixfile(myfilestr):
    infile = open(myfilestr, "rb" )
    instr = infile.read()
    infile.close()
    outstr = instr.replace( "\r\n", "\n" ).replace( "\r", "\n" ).replace( "\n", "\r\n" )

    if len(outstr) == len(instr):
        print 'same length'
        return 

    outfile = open(myfilestr, "wb" )
    outfile.write( outstr )
    outfile.close()


def readsansfile(myfilestr):
    fixfile(myfilestr)
    myfile=open(myfilestr,"r")
    count=0
    header=[]
    myFlag=True
    data=[]
    while myFlag:
        try:
            returnline=['']
            if count <17:
                currline=get_tokenized_line(myfile,returnline)
                header.append(returnline)
                print len(currline), currline
                if count==0:
                    mytimestampstr=currline[1]+' '+currline[2]
                    mydatetime=datetime.strptime(mytimestampstr[1:-2], "%d-%b-%Y %H:%M:%S")
                    mytimestamp=time.mktime(mydatetime.timetuple())
                if count==2:
                    monitor=float(currline[1])          
            else:
                currline=get_tokenized_line(myfile,returnline,splitchar=',')
                print len(currline), currline
                if currline in ['',['']]:
                    myFlag=False
                else:
                    #data.append(np.array(currline[:-1],'Float64'))
                    data=np.concatenate((data,np.array(currline[:-1],'Float64')))
            count=count+1
        except:
            print 'done'
            myFlag=False        
    res={}
    res['header']=header
    res['monitor']=monitor
    res['timestamp']=mytimestamp
    res['data']=np.array(data)
    return res
    #Epoch=float(tokenized[1])
    #timeobj=mx.DateTime.DateTimeFromTicks(ticks=Epoch) #what I originally used
    #timeobj=datetime.fromtimestamp(Epoch)
    

def readfiles(ppfile):
    myfiledir=r"c:\bfosans\data\grasp"
    myfiledir=r"c:\bfosans\data"
    #ppfile=150
    myend='.GSP'
    myname="BFORV"+str(ppfile)
    myfilestr=os.path.join(myfiledir,myname)
    myfilestr=myfilestr+myend
    #return myfilestr
    result=readsansfile(myfilestr)
    print result
    return result
    
def read_driver(myfilestr):
    fixfile(myfilestr)
    myfile=open(myfilestr,"r")
    count=0
    myFlag=True
    res=[]
    #++    -+    +-    --    
    #* denotes absent scan
    results=[]
    while myFlag:
        try:
            currline=get_tokenized_line(myfile)
            resdict={}
            if currline[0]!='*':
                resdict['off_off']=currline[0]
            if currline[1]!='*':
                resdict['on_off']=currline[1]
            if currline[2]!='*':
                resdict['off_on']=currline[2]
            if currline[3]!='*':
                resdict['off_off']=currline[3]
            results.append(resdict)     
        except:
            print 'done'
            myFlag=False
        
    return results
            
        
        
        
    
    

if __name__=="__main__":
    myfiledir=r"c:\bfosans"
    myname="Pnki_BFORV_Scans.txt"
    myname="Ppki_BFORW_Scans.txt" # big endian??
    myname="p2.txt"
    myfilestr=os.path.join(myfiledir,myname)
    results=read_driver(myfilestr)
    keys=['off_off','off_on','on_off','on_on']
    i=0
    res={}
    for key in keys:
        try:
            fp=results[i][key]
            res[key]=readfiles(fp)
            print 'read'
        except:
            print key, 'does not exist for',i
    print res
