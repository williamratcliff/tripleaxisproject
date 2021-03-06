import numpy as N
import pylab
from utilities import  scriptutil as SU
import re,os,sys
from utilities import readncnr3 as readncnr
import copy
from reflectometry.reduction.rebin import rebin
pi=N.pi



#class 

def cshift(l, offset):
   offset %= len(l)
   return N.concatenate((l[-offset:], l[:-offset]))



class Psd:
    def __init__(self,left=0,center=21,right=37):
        self.left=left
        self.right=right
        self.center=center
        return

class Stitch:
    def __init__(self,data,ch_a4_str,ch_eff_str,psd,outputwidth=0.31,masked=[]):
        
        ch_a4=N.loadtxt(ch_a4_str, unpack=True)
        ch_a4=ch_a4.T.flatten()
        ch_eff=N.loadtxt(ch_eff_str, unpack=True)
        ch_eff=ch_eff.T.flatten()
        

        #ch_eff=N.ones(ch_eff.shape)
        self.ch_eff_orig=ch_eff
        ch_eff[masked]=0.0
        self.ch_eff=ch_eff
        self.data=data
        self.psd=psd
        ch_space=ch_a4[psd.center]-ch_a4
        print 'ch_space1', ch_space
        ch_space=N.concatenate((ch_space,[ch_space[-1]]))
        print 'ch_space2',ch_space
        print 'shifted',cshift(ch_space,1)
        ch_boundary=ch_space+0.5*(-ch_space+cshift(ch_space,1))
        ch_boundary[0]=ch_space[0]+0.5*(ch_space[0]-ch_space[1])
        ch_boundary[-1]=ch_space[-2]-0.5*(ch_space[-3]-ch_space[-2])
        self.ch_space=ch_space
        self.ch_boundary=ch_boundary

        
        self.outputwidth=outputwidth
        self.a4=N.array(data.data['a4'],'float64')
        self.a4_begin=self.a4[0]-abs(ch_space[psd.right]) #right is lower a4
        self.a4_end=self.a4[-1]+abs(ch_space[psd.left]) #left is higher a4
        #self.a4_begin=30.0
        self.output_npts=int(N.round(N.absolute(self.a4_end-self.a4_begin)/self.outputwidth))
        #self.output_a4=N.arange(N.min([self.a4_begin,self.a4_end]),N.max([self.a4_begin,self.a4_end])+outputwidth,outputwidth)
        self.output_a4=N.linspace(N.min([self.a4_begin,self.a4_end]),N.max([self.a4_begin,self.a4_end]),self.output_npts+1)
        #print 'output',self.output_a4
        #self.output_npts=self.output_a4.shape[0]
        #print self.output_npts
        #print self.output_a4.shape
        print 'begin',self.a4_begin
        print 'end',self.a4_end 
        print 'boundaries',ch_boundary

        if 0:
            pylab.plot(ch_boundary,'s')
            #pylab.plot(ch_a4,'s')
            #pylab.plot(ch_eff,'s')
            pylab.show()
            sys.exit()

        
        #print 'a4',self.a4       
        detectors=self.data.metadata['count_info']['AnalyzerDetectorDevicesOfInterest'.lower()]
        self.corrected=False
        self.counts_ending=''
        self.errors_ending=''
        self.detectors=detectors
        if data.data.has_key(detectors[0]+'_corrected'):
            self.corrected=True
            self.counts_ending='_corrected'
            self.errors_ending='_errs_corrected'
        i=0
        z_in=[]
        dz_in=[]
        for detector in detectors:
            z_in.append(N.array(self.data.data[detector+self.counts_ending])*ch_eff[i])
            if self.errors_ending!='':
                dz_in.append(N.array(self.data.data[detector+self.errors_ending])*ch_eff[i])
            else:
                dz_in.append(N.sqrt(N.array(self.data.data[detector+self.counts_ending]))*ch_eff[i])
            #print detector+self.counts_ending,N.array(self.data.data[detector+self.counts_ending])
            i=i+1
        self.data_eff=N.array(z_in)#.T
        self.data_err_eff=N.array(dz_in)#.T
        #print self.data_eff
        #print self.data_err_eff
        print 'shape',self.data_eff.shape
        #for i in range(13):
        #    pylab.plot(self.data_eff[:,i],'s')
        
        if 0:
            ch_space2=ch_a4[psd.center]-ch_a4
            print ch_space2
            print self.a4
            print masked
            a4size=self.a4.shape[0]
            for i,a in enumerate(self.a4):
                #pylab.plot(a+ch_space2,self.data_eff[:,i],'s')
                if i in [a4size-3,a4size-2,a4size-1]:
                    pylab.plot(a+ch_space2,self.data_eff[:,i],'s')
            pylab.show()
            sys.exit()
        

        
        
    def stitch(self):
        data=self.data
        #data_plus_count=N.zeros(self.output_a4.shape,'float64')
        mon_in=N.ones(self.ch_eff.shape,'float64')
        i=0
        ch_eff=self.ch_eff
        #ch_eff[[3,7,8]]=0
        for detector in self.detectors:
            if ch_eff[i] < 1e-2:
                mon_in[i]=0
            i=i+1
        #print mon_in
        #output_tmp=N.zeros(self.output_a4.shape,'float64')
        #dz_output_tmp=N.zeros(self.output_a4.shape,'float64')
        #output_mon=N.zeros(self.output_a4.shape,'float64')
        output_data=0.0*self.output_a4[0:-1]
        output_data_err=N.zeros(output_data.shape,'float64')
        data_norm=N.zeros(output_data.shape,'float64')
        mon_in=z_in=N.flipud(mon_in)
        for i in range(self.a4.shape[0]):
            #ch_boundary=self.ch_boundary+self.a4[i]
            z_in = self.data_eff[:,i]#*ch_eff
            dz_in = self.data_err_eff[:,i]#*ch_eff
            #print z_in
            #print dz_in
            ch_boundary=N.flipud(self.ch_boundary+self.a4[i])
            print 'inline boundary', ch_boundary
            z_in=N.flipud(z_in)
            dz_in=N.flipud(dz_in)
            #;output_data_left=reverse(output_data_left)
            output_tmp=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],z_in[self.psd.left:self.psd.right],self.output_a4))
            dz_output_tmp=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],dz_in[self.psd.left:self.psd.right],self.output_a4))
            output_mon=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],mon_in[self.psd.left:self.psd.right],self.output_a4))
            #drebin,ch_left[range[0]:range[2]+1],z_in[range[0]:range[2]],dz_in[range[0]:range[2]],output_data_left,output_tmp,dz_output_tmp,/histogram,/to_histogram,err=err,emsg=emsg
            
            #print ch_boundary[self.psd.left:self.psd.right+1].tolist()
            #print self.output_a4
            #print z_in[self.psd.left:self.psd.right].tolist()
            #print self.output_a4.tolist()
            
            #print mon_in
            #print 'output_mon',output_mon
            #print, emsg
            #drebin,ch_left[range[0]:range[2]+1],mon_in[range[0]:range[2]],mon_in[range[0]:range[2]],output_data_left,output_mon,dz_output_mon,/histogram,/to_histogram,err=err,emsg=emsg
            print output_data.shape
            print output_tmp.shape
            #print output_mon.shape
            #print dz_output_tmp.shape
            #print output_tmp^2
            output_data=output_data+output_tmp#*output_mon
            output_data_err=output_data_err+dz_output_tmp**2#*output_mon**2
            #print 'output_mon',output_mon
            data_norm=data_norm+output_mon
            #print 'data_norm',data_norm
            
            
            
        idx=N.where(data_norm>0)
        self.output_data=output_data[idx]/data_norm[idx]
        self.output_data_err=N.sqrt(output_data_err[idx])/data_norm[idx]
        self.output_a4=self.output_a4[idx]
        #print 'output_a4', self.output_a4
        print 'data_norm', data_norm
        #print 'output',self.output_data
        #print 'norm'
        #print data_norm
        #print data_norm.shape
        #print output_mon.shape
        #print dz_output_tmp.shape
        #print output_data.shape
        #print output_tmp.shape
        
    def writefile(self,myoutfilestr):
        myoutfile=open(myoutfilestr,'wt')
        a4=self.output_a4[0:-1]
        counts=self.output_data
        for i in range(a4.shape[0]):
            myoutfile.write('%3.3f %3.3f\n'%(a4[i],counts[i]))
        myoutfile.close()

if __name__=='__main__':
        mydirectory=r'C:\13165\13165\data'
        myfilebase='PSD_A4_SpacingApr1108.dat'
        myend='bt7'
        #myfilestr=r'C:\13165\PSD_A4_SpacingApr1108.dat'
        #myfilestr=r'C:\13165\PSD_Channel_April1208.dat'
        ch_a4_str=r'C:\13188\PSD_A4_Spacing_Jun2008_Ana_85mm_80minRC.dat'        
        ch_eff_str=r'c:\13188\PSD_Channel_Eff_Jun202008_Ana_85mm_80minRC.dat'
        ch_eff_str=r'c:\13188\eff.dat'
        
        #jun
        ch_a4_str=r'C:\jun\PSD_A4_SpacingApr1108.dat'        
        ch_eff_str=r'c:\jun\PSD_Channel_EffJune0908.dat'
        
        
        #clarina LiFeas
        ch_a4_str=r'C:\NaFeAs2\PSD_spacing_40minradial_Feb142009.dat'
        ch_eff_str=r'C:\NaFeAs2\PSD_efficiency_40minradial_Feb142009.dat'

        mypsd=Psd(center=23)

        
        #me
        #mydirectory2=r'C:\13176\data'
        #myfilestr=os.path.join(mydirectory2,'CeOFeAs57256.bt7.out')
        #ying
        mydirectory2=r'C:\13188'
        myfilestr=os.path.join(mydirectory2,'NdOFeAs58081.bt7.out')
        myfilestr=os.path.join(mydirectory2,'NdOFeAs58075.bt7')
        myfilestr=os.path.join(mydirectory2,'NdOFeAs58054.bt7')
        #myfilestr=os.path.join(mydirectory2,'plastic_80RC58038.bt7')
        #myfilestr=os.path.join(mydirectory2,'plastic57943.bt7')
        
        #jun
        mydirectory2=r'C:\jun'
        #myfilestr=os.path.join(mydirectory2,'LiFeA002.bt9')
        myfilestr=os.path.join(mydirectory2,'CeOFeAs859493.bt7')
        
        #clarina LiFeAs
        mydirectory2=r'C:\NaFeAs2'
        myfilestr=os.path.join(mydirectory2,'NaFeAs68019.bt7')
        myfilestr=os.path.join(mydirectory2,'plastic69167.bt7')
        
        
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        mystitcher=Stitch(mydata,ch_a4_str,ch_eff_str,mypsd,masked=[])        
        mystitcher.stitch()
        
        #print mystitcher.data_eff.shape
        #print mystitcher.data_eff[0,:]
        #print mystitcher.ch_eff
        #print mystitcher.detectors
        #print mystitcher.data.data[mystitcher.detectors[-1]+'_corrected']
        myoutstr=myfilestr+'.stitched'
        mystitcher.writefile(myoutstr)
        tdata=mystitcher.output_data
        if 1:
            #print 'data_eff', mystitcher.data_eff
            #pylab.plot(mystitcher.data_eff,'s')
            #pylab.errorbar(mystitcher.output_a4[0:-1],mystitcher.output_data,mystitcher.output_data_err,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
            pylab.errorbar(mystitcher.output_a4,mystitcher.output_data,mystitcher.output_data_err,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
            #print 'a',mystitcher.output_a4[0:-1]
        if 0:
            #myfilestr=os.path.join(mydirectory2,'CeOFeAs57257.bt7.out')
            #myfilestr=os.path.join(mydirectory2,'NdOFeAs58082.bt7.out')
            #myfilestr=os.path.join(mydirectory2,'NdOFeAs58076.bt7.out')
            myfilestr=os.path.join(mydirectory2,'NdOFeAs58047.bt7.out')
            
            mydatareader=readncnr.datareader()
            mydata=mydatareader.readbuffer(myfilestr)
            mystitcher=Stitch(mydata,ch_a4_str,ch_eff_str,mypsd,masked=[6])
            mystitcher.stitch()
            myoutstr=myfilestr+'.stitched'
            mystitcher.writefile(myoutstr)
            print tdata.shape
            print mystitcher.output_data.shape
            #print 'b',mystitcher.output_a4[0:-1]
            ddata=mystitcher.output_data#-tdata
            pylab.errorbar(mystitcher.output_a4[0:-1],mystitcher.output_data,mystitcher.output_data_err,marker='s',linestyle='None',mfc='red',mec='black',ecolor='black')
            
            #pylab.axis([33,36,])
    
            
        
        if 1:
            #print mystitcher.output_a4[0:-1]
            print mystitcher.output_a4[0]
            print mystitcher.output_a4[-1]
            pylab.show()
        
        
        