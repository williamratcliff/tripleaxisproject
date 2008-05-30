import numpy as N
import pylab
from utilities import  scriptutil as SU
import re,os,sys
from utilities import readncnr3 as readncnr
import copy
from reflectometry.reduction import rebin
pi=N.pi



#class 

def cshift(l, offset):
   offset %= len(l)
   return N.concatenate((l[-offset:], l[:-offset]))

def rebin(x,I,xo,Io=None):
    """
    Rebin a vector.
    
    x are the existing bin edges
    xo are the new bin edges
    I are the existing counts (one fewer than edges)
    
    Io will be used if present, but be sure that it is a contiguous
    array of the correct shape and size.
    """
    x,I,xo = [_input(v) for v in (x,I,xo)]
    Io = _output(Io,[len(xo)-1])
    _reduction.rebin(x,I,xo,Io)
    return Io

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
        self.output_a4=N.arange(N.min([self.a4_begin,self.a4_end]),N.max([self.a4_begin,self.a4_end]),outputwidth)
        
        print self.output_npts
        print self.output_a4
        print self.a4_end        
        
    def stitch(self):
        data=self.data
        for i in range(self.a4.shape[0]):
            ch_boundary=self.ch_boundary+self.a4[i]
            
        
        


if __name__=='__main__':
        mydirectory=r'C:\13165\13165\data'
        myfilebase='PSD_A4_SpacingApr1108.dat'
        myend='bt7'
        flist=[]
        rlist=[56421,56420,56419]
        #myfilestr=r'C:\13165\PSD_A4_SpacingApr1108.dat'
        myfilestr=r'C:\13165\PSD_A4_Spacing.dat'
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
        
        #print ch_space
        #print ch_boundary
        myfilestr=os.path.join(mydirectory,'CeOFeAs56685.bt7')
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        mystitcher=Stitch(mydata,ch_space,ch_boundary,ch_eff,mypsd)
        outdata=mystitcher.stitch()
        
        
        
        
#        pylab.plot(N.arange(ch_eff.size),ch_eff,'ro')
#        pylab.show()