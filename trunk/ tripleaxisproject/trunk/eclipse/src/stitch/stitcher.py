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
    def __init__(self,left=0,center=21,right=48):
        self.left=left
        self.right=right
        self.center=center
        return

class Stitch:
    def __init__(self,data,ch_space,ch_boundary,ch_eff,psd,outputwidth=0.1):
        self.ch_boundary=ch_boundary
        self.ch_eff=ch_eff
        self.ch_space=ch_space
        self.data=data
        self.psd=psd
        self.outputwidth=outputwidth
        self.a4=N.array(data.data['a4'],'float64')
        self.a4_begin=self.a4[0]-ch_space[psd.left]
        self.a4_end=self.a4[-1]-ch_space[psd.right]
        self.output_npts=int(N.round(N.absolute(self.a4_end-self.a4_begin)/self.outputwidth))
        self.output_a4=N.arange(N.min([self.a4_begin,self.a4_end]),N.max([self.a4_begin,self.a4_end])+outputwidth,outputwidth)
        
        #print self.output_npts
        #print self.output_a4.shape
        #print self.a4_begin
        #print self.a4_end        
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
        
        #print self.data_eff[:,0]
        #print ch_eff
        #print ch_eff[0]
        #print z_in^2
        
        #print self.data_eff.shape
        #print z_in.shape
        
        
    def stitch(self):
        data=self.data
        #data_plus_count=N.zeros(self.output_a4.shape,'float64')
        mon_in=N.ones(self.ch_eff.shape,'float64')
        i=0
        ch_eff=self.ch_eff
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
        data_norm=N.ones(output_data.shape,'float64')
        for i in range(self.a4.shape[0]):
            ch_boundary=self.ch_boundary+self.a4[i]
            z_in = self.data_eff[:,i]#*ch_eff
            dz_in = self.data_err_eff[:,i]#*ch_eff
            #print z_in
            #print dz_in
            ch_boundary=N.flipud(ch_boundary)
            z_in=N.flipud(z_in)
            dz_in=N.flipud(dz_in)
            #;output_data_left=reverse(output_data_left)
            output_tmp=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],z_in[self.psd.left:self.psd.right],self.output_a4))
            dz_output_tmp=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],dz_in[self.psd.left:self.psd.right],self.output_a4))
            output_mon=N.array(rebin(ch_boundary[self.psd.left:self.psd.right+1],mon_in[self.psd.left:self.psd.right],self.output_a4))
            #drebin,ch_left[range[0]:range[2]+1],z_in[range[0]:range[2]],dz_in[range[0]:range[2]],output_data_left,output_tmp,dz_output_tmp,/histogram,/to_histogram,err=err,emsg=emsg
            
            print ch_boundary[self.psd.left:self.psd.right+1].tolist()
            print self.output_a4
            #print z_in[self.psd.left:self.psd.right].tolist()
            print self.output_a4.tolist()
            
            #print mon_in
            print output_mon
            #print, emsg
            #drebin,ch_left[range[0]:range[2]+1],mon_in[range[0]:range[2]],mon_in[range[0]:range[2]],output_data_left,output_mon,dz_output_mon,/histogram,/to_histogram,err=err,emsg=emsg
            #print output_data.shape
            #print output_tmp#.shape
            #print output_mon.shape
            #print dz_output_tmp.shape
            #print output_tmp^2
            output_data=output_data+output_tmp*output_mon
            output_data_err=output_data_err+dz_output_tmp**2*output_mon**2
            data_norm=data_norm+output_mon
        self.output_data=output_data/data_norm
        self.output_data_err=N.sqrt(output_data_err)/data_norm
        print 'norm'
        print data_norm
        #print data_norm.shape
        #print output_mon.shape
        #print dz_output_tmp.shape
        #print output_data.shape
        #print output_tmp.shape
        


if __name__=='__main__':
        mydirectory=r'C:\13165\13165\data'
        myfilebase='PSD_A4_SpacingApr1108.dat'
        myend='bt7'
        flist=[]
        rlist=[56421,56420,56419]
        myfilestr=r'C:\13165\PSD_A4_SpacingApr1108.dat'
        #myfilestr=r'C:\13165\PSD_A4_Spacing.dat'
        ch_a4=N.loadtxt(myfilestr, unpack=True)
        ch_a4=ch_a4.T.flatten()
        myfilestr=r'C:\13165\PSD_Channel_April1208.dat'
        ch_eff=N.loadtxt(myfilestr, unpack=True)
        ch_eff=ch_eff.T.flatten()
        #ch_space+0.5*(-ch_space+shift(ch_space,1))+var(var_pos,i)
        #print ch_eff
        mypsd=Psd(center=23)
        #print ch_a4
        #print (ch_a4-ch_a4[mypsd.center])
        ch_space=ch_a4[mypsd.center]-ch_a4
        ch_space=N.concatenate((ch_space,[ch_space[-1]]))
        #ch_boundary=(ch_a4-ch_a4[mypsd.center])+0.5*N.concatenate((N.array([0]),N.diff(ch_a4)))
        #ch_boundary[0]=(ch_a4[0]-ch_a4[mypsd.center])-0.5*(ch_a4[1]-ch_a4[0])
        ch_boundary=ch_space+0.5*(-ch_space+cshift(ch_space,1))
        ch_boundary[0]=ch_space[0]+0.5*(ch_space[0]-ch_space[1])
        ch_boundary[-1]=ch_space[-2]-(ch_space[-3]-ch_space[-2])
        
        #print N.array(ch_space[-1])
        
        #print ch_a4
        #print ch_space
        #print ch_boundary
        myfilestr=os.path.join(mydirectory,'CeOFeAs56685.bt7')
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        mystitcher=Stitch(mydata,ch_space,ch_boundary,ch_eff,mypsd)
        mystitcher.stitch()
        print mystitcher.output_a4
        #print mystitcher.output_data.shape
        #pylab.plot(mystitcher.output_a4[0:-1],mystitcher.output_data)
        pylab.errorbar(mystitcher.output_a4[0:-1],mystitcher.output_data,mystitcher.output_data_err,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
        pylab.show()
        
        
        
        
#        pylab.plot(N.arange(ch_eff.size),ch_eff,'ro')
#        pylab.show()