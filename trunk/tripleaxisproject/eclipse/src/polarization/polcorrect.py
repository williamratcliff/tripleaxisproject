import numpy as N
#import pylab
import math
import readncnr2 as readncnr
import writebt7
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
#cellP="PcellBT7Jan072008.txt"
#cellA="AcellBT7Jan72008.txt"
#cellP="PcellBT7Jan92008.txt"
#cellA="AcellBT7Jan92008.txt"


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

mypolcorrect = N.ctypeslib.load_library('polarization2.dll', '.')




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

def calc_energy(angle,dspace):
    anglerad=myradians(angle)
    tau=2*N.pi/dspace
    k=tau/N.sqrt(2-2*N.cos(anglerad))
    energy=2.072142*k*k
    return energy

class polarization_correct:
    def __init__(self,files,cell):
        self.mydata={}
        self.counts={}
        self.errors={}
        ei=[]
        ef=[]
        self.timestamp={}
        self.files=files
        self.cell=cell
        for key,myfilestr in files.iteritems():
            mydatareader=readncnr.datareader()
            self.mydata[key]=mydatareader.readbuffer(myfilestr)
            #print mydata[key].metadata
            self.counts[key]=N.array(self.mydata[key].data[self.mydata[key].metadata['count_info']['signal']])
            self.errors[key]=N.sqrt(self.counts[key])
            self.timestamp[key]=N.array(self.mydata[key].data['timestamp'])
            #TODO currently, we assume that the files are taken at the same points--is this safe?
            self.length=self.counts[key].shape[0]
            a2=N.array(self.mydata[key].data['a2'])
            if self.mydata[key].metadata['count_info']['analyzerdetectormode'].lower()=='diffdet':
                #a6=N.array(mydata[key].data['A5'])*2.0
                a6=a2
            else:
                a6=N.array(self.mydata[key].data['a6'])
            #TODO who's bright idea was it to have A6 listed as "IN"
            #TODO we also assume that the energies will be the same!!!
            dmono=self.mydata[key].metadata['dspacing']['monochromator_dspacing']
            dana=self.mydata[key].metadata['dspacing']['analyzer_dspacing']
            #print a6
            ei.append(calc_energy(a2,dmono))
            ef.append(calc_energy(a6,dana))
            #print self.length
        self.ei=N.array(ei,'float64')
        self.ef=N.array(ef,'float64')
        #print mydata[key].data.keys()
        return
    def output(self,outputfile=None):
        keys=['pp','mm','pm','mp']
        s=''
        if outputfile!=None:
            f=open(outputfile,'wt')
        for i in range(self.length):
            s=s+'%f %f '%(self.ei[0][i],self.ef[0][i])
            for key in keys:
                if self.counts.has_key(key):
                    s=s+'%f '%(self.counts[key][i],)
                    s=s+'%f '%(self.timestamp[key][i],)
                else:
                    s=s+'* * '
            s=s+'\n'
            print i
        if outputfile==None:
            print s
        else:
            f.write(s)
        if outputfile!=None:
            f.close()
        return

    def savefiles(self):
        mywriter=writebt7.datawriter()
        for key,myfilestr in self.files.iteritems():
            mywriter.write(myoutfilestr=myfilestr+'.out',mydata=self.outdata[key]) 
        return

    
    def correct(self):
        keys=['pp','mm','pm','mp']
        pbinput=PBindata()
        pboutput=PBoutdata()
        pbinput.Ei=self.ei.ctypes.data_as(c_double_p)
        pbinput.Ef=self.ef.ctypes.data_as(c_double_p)
        for key in keys:
            if self.counts.has_key(key):
                if key=='pp':
                    pbinput.tpp=self.timestamp[key].astype('uint32').ctypes.data_as(c_ulong_p)
                    pbinput.Cpp=self.counts[key].ctypes.data_as(c_double_p)
                    pbinput.Epp=self.errors[key].ctypes.data_as(c_double_p)
                if key=='mm':
                    pbinput.tmm=self.timestamp[key].astype('uint32').ctypes.data_as(c_ulong_p)
                    pbinput.Cmm=self.counts[key].ctypes.data_as(c_double_p)
                    pbinput.Emm=self.errors[key].data_as(c_double_p)
                if key=='pm':
                    mytemp_pm=self.timestamp[key].astype('uint32')
                    pbinput.tpm=mytemp_pm.ctypes.data_as(c_ulong_p)
                    #print 't0 ',mytemp_pm[0]
                    pbinput.Cpm=self.counts[key].ctypes.data_as(c_double_p)
                    pbinput.Epm=self.errors[key].ctypes.data_as(c_double_p)
                    #print 'shape ',self.counts[key].shape
                if key=='mp':
                    mytemp_mp=self.timestamp[key].astype('uint32')
                    pbinput.tmp=mytemp_mp.ctypes.data_as(c_ulong_p)
                    pbinput.Cmp=self.counts[key].ctypes.data_as(c_double_p)
                    pbinput.Emp=self.errors[key].ctypes.data_as(c_double_p)
##                tkey='t'+key
##                pbinput[tkey]=
##                ckey='C'+key
##                pbinput[ckey]=self.counts[key].data
##                ekey='E'+key
##                pbinput[ekey]=(N.sqrt(self.counts[key])).data
            else:
##                tkey='t'+key
##                pbinput[tkey]=N.array([],'float64').ctypes.data_as(c_ulong_p)
##                ckey='C'+key
##                pbinput[ckey]=N.array([],'float64').data
##                ekey='E'+key
##                pbinput[ekey]=N.array([],'float64').data
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
               
            
        #print 'mylength ',self.length
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
        pbflags=PBflags()

        fMonitorCorrect=0
        fPolMonitorCorrect=0
        if self.mydata[key].metadata['count_info']['count_type']=='monitor':
            fPolMonitorCorrect=1
            #print 'monitor'
        fDebug=0
        fSimFlux=0
        fSimDeviate=0
        #-- ++ +- -+ #TODO Check with Ross on the order


        pbflags.MonitorCorrect=fMonitorCorrect
        pbflags.PolMonitorCorrect=fPolMonitorCorrect
        pbflags.Debug=fDebug
        pbflags.SimFlux=fSimFlux
        pbflags.SimDeviate=fSimDeviate
        pbflags.CountsEnable=fCountsEnable.ctypes.data_as(c_int_p)
        pbflags.CountsAdd1=fCountsAdd1.ctypes.data_as(c_int_p)
        pbflags.CountsAdd2=fCountsAdd2.ctypes.data_as(c_int_p)
        pbflags.Sconstrain=fSconstrain.ctypes.data_as(c_int_p)
        pbflags.Spp=fSpp.ctypes.data_as(c_double_p)
        pbflags.Smm=fSmm.ctypes.data_as(c_double_p)
        pbflags.Spm=fSpm.ctypes.data_as(c_double_p)
        pbflags.Smp=fSmp.ctypes.data_as(c_double_p)
        mypolcorrect.PBcorrectData(cellP,cellA,ctypes.byref(pbflags),self.length,ctypes.byref(pbinput),ctypes.byref(pboutput))
 
        pbflags.MonitorCorrect=0#fMonitorCorrect
        pbflags.PolMonitorCorrect=fPolMonitorCorrect
        #print 'count type ',self.mydata[key].metadata['count_info']['count_type']
        #print fPolMonitorCorrect
        pbflags.MonoSelect=1
        pbflags.Debug=0#fDebug
        pbflags.SimFlux=0#fSimFlux
        pbflags.SimDeviate=0#fSimDeviate
        pbflags.NoNegativeCS=0
        pbflags.HalfPolarized=0
        pbflags.CountsEnable=[0,0,1,1]
        pbflags.CountsAdd1=[0,0,0,0]
        pbflags.CountsAdd2=[0,0,0,0]
        pbflags.Sconstrain=[1,1,0,0]
        pbflags.Spp=[0,0,0,0]
        pbflags.Smm=[0,0,0,0]
        bob=N.empty((1,3),'int32')
        pbflags.Spm=[0,0,0,0]
        pbflags.Smp=[0,0,0,0]
   
        mypolcorrect.PBcorrectData(self.cell,pbflags._pointer,self.length,ctypes.byref(pbinput),ctypes.byref(pboutput))
#int PBcorrectData(char *PCellFile, char *ACellFile, PBflags *flgs,
#		  int npts, PBindata *in, PBoutdata *out) ;
        #print pboutput.Spm
        #print Smp[0]
        #print Spm[0]
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
        
        #append corrected counts to our dataset
        self.outdata={}
        for key in keys:
            if self.counts.has_key(key):
                self.outdata[key]=copy.deepcopy(self.mydata[key])
                #self.mydata[key].data[self.mydata[key].metadata['signal']]
                newfield=self.mydata[key].metadata['count_info']['signal']+'_corrected'
                self.outdata[key].data[newfield]=corrected_counts['S'+key]
                newfield=self.mydata[key].metadata['count_info']['signal']+'_errs_corrected'
                self.outdata[key].data[newfield]=corrected_counts['E'+key]
            
        #print Epm[0]
        #print Epp[0]
        return corrected_counts


def readscript(myfilestr):
    myfile = open(myfilestr, 'r')
    mycount=0
    #myfile.close()
    #exit()
    #myfile = open(myfilestr, 'r')
    returnline=['']
    catalog=[]
    cellfiles=[]
    absolute=True
    inblock=False
    while 1:
        lineStr = myfile.readline()
        if not(lineStr):
            break
        strippedLine=lineStr.rstrip().lower()
        tokenized=strippedLine.split()
        #print 'tokenized ', tokenized
        if tokenized[0]==[]:
            pass
        elif tokenized[0]=='#absolute':
            #print 'absolute'
            absolute=True
            mydirectory=os.curdir
        elif tokenized[0]=='#directory':
            #print 'directory'
            if os.path.isdir(mydirectory):
                mydirectory=tokenized[1]
            else:
                print '%s is not an existing directory!!!'%(tokenized[1],)
                break
        elif tokenized[0]=='#begin':
            inblock=True
            #print 'begin '
            files={}
        elif tokenized[0]=='#end':
            #print 'end'
            catalog.append(files)
            print 'correcting %s using cellfile %s'%(files,cellfile)
            mypolcor=polarization_correct(files,cellfile)
            corrected_counts=mypolcor.correct()
            mypolcor.savefiles()
            print 'corrected'
            if 0:
                key='pm'
                #pylab.subplot(2,2,1+mycount)
                #pylab.title(key)
                mydatac=mypolcor.mydata
                #pylab.errorbar(mydatac[key].data['qx'],mydatac[key].data['detector'],N.sqrt(mydatac[key].data['detector']),
                #    marker='s',mfc='blue',linestyle='None')
                #pylab.errorbar(mydatac[key].data['qx'],corrected_counts['Spm'],corrected_counts['Epm'], marker='s',mfc='red',linestyle='None')
                #print 'pm'
                #print corrected_counts['Spm']
                key='mp'
                #pylab.subplot(2,2,2+mycount)
                #pylab.title(key)
                #print 'mp'
                #print corrected_counts['Smp']
                #pylab.errorbar(mydatac[key].data['qx'],mydatac[key].data['detector'],N.sqrt(mydatac[key].data['detector']),
                #    marker='s',mfc='blue',linestyle='None')
                #pylab.errorbar(mydatac[key].data['qx'],corrected_counts['Smp'],corrected_counts['Emp'], marker='s',mfc='red',linestyle='None')
                mycount=mycount+2
            inblock=False
        elif inblock==True:
            #print 'inblock'
            toksplit=tokenized[0].split('=')
            #check to make sure that there actually is a file specified!
            if len(toksplit)==2:
                filetok=os.path.join(mydirectory,toksplit[1])
                if os.path.isfile(filetok):
                    files[toksplit[0]]=filetok
                else:
                    print filetok+' does not exist!!!'
                    break
            else:
                pass
        elif tokenized[0]=='#cell'.lower():
            #print 'acellfile'
            cellfile=os.path.join(mydirectory,tokenized[1])
            if os.path.isfile(cellfile):
                cellfiles.append(cellfile)
            else:
                print '%s does not exist!'%(cellfile,)
                break
        #elif tokenized[0]=='#Pcell'.lower():
        #    #print 'pcellfile'
        #    pcellfile=os.path.join(mydirectory,tokenized[1])
        #    if os.path.isfile(pcellfile):
        #        pcellfiles.append(pcellfile)
        #    else:
        #        print '%s does not exist!'%(acellfile,)
        #        break
    print 'done'
    myfile.close()
    print 'closed'

if __name__=="__main__":
#    myfilestr_on_off=r'c:\bifeo3xtal\jan8_2008\9175\fieldscansplusminus53566.bt7'
#    myfilestr_off_on=r'c:\bifeo3xtal\jan8_2008\9175\fieldscanminusplus53567.bt7'
    myfilestr_on_off=r'c:\bifeo3xtal\jan8_2008\9175\fieldscansplusminusreset53630.bt7'
    myfilestr_off_on=r'c:\bifeo3xtal\jan8_2008\9175\fieldscanminusplusreset53631.bt7'
    myscriptstr=r'c:\tripleaxisproject2\polcor\bifeo3.polcor'

    if len(sys.argv)!=2:
        print 'the usage of this program is polcorrect scriptfile'
        print 'scriptfiles consist of directives set with a #'
        print '#absolute on a line tells me that you will use absolute paths for files'
        print 'if you do not specify a path, it will default to the current working directory'
        print '#directory directory_name tells me that you will use directory_name for your files'
        print '#cell cell_file name gives me the name of the cell file'
        print '#begin tells me that you are about to give a list of files that should be grouped together for reduction'
        print 'within this block, you can use pm=filename on a single line'
        print 'the valid channels are pp,pm,mp,mm  with polarizer in being denoted by p and out denoted by m--given in stream order'
        print '#end tells me that the block is finished.  Do not forget this directive!'
        exit()
    myscriptstr=sys.argv[1]
    if os.path.isfile(myscriptstr):
        readscript(myscriptstr)
    else:
        print 'The file you called me with (%s) does not exist'%(myscriptstr,)
    #pylab.show()
    exit()
