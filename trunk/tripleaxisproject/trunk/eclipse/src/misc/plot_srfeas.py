import numpy as N
import utilities.readncnr3 as readncnr
import sys, os
import pylab
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import NullFormatter, MultipleLocator

from math import fmod
import numpy

import matplotlib.cbook as cbook
import matplotlib.transforms as transforms
import matplotlib.artist as artist
import matplotlib.patches as patches
from matplotlib.path import Path


class Ring(patches.Patch):
    """
    Ring patch.
    """
    def __str__(self):
        return "Ring(%g,%g,%g,%g)"%(self.r1,self.r2,self.theta1,self.theta2)

    def __init__(self,
                 center=(0,0),
                 r1=0,
                 r2=None,
                 theta1=0,
                 theta2=360,
                 **kwargs
                 ):
        """
        Draw a ring centered at *x*, *y* center with inner radius *r1* and
        outer radius *r2* that sweeps *theta1* to *theta2* (in degrees).

        Valid kwargs are:

        %(Patch)s
        """
        patches.Patch.__init__(self, **kwargs)
        self.center = center
        self.r1, self.r2 = r1,r2
        self.theta1, self.theta2 = theta1,theta2

        # Inner and outer rings are connected unless the annulus is complete
        delta=abs(theta2-theta1)
        if fmod(delta,360)<=1e-12*delta:
            theta1,theta2 = 0,360
            connector = Path.MOVETO
        else:
            connector = Path.LINETO

        # Form the outer ring
        arc = Path.arc(theta1,theta2)

        if r1 > 0:
            # Partial annulus needs to draw the outter ring
            # followed by a reversed and scaled inner ring
            v1 = arc.vertices
            v2 = arc.vertices[::-1]*float(r1)/r2
            v = numpy.vstack([v1,v2,v1[0,:],(0,0)])
            c = numpy.hstack([arc.codes,arc.codes,connector,Path.CLOSEPOLY])
            c[len(arc.codes)]=connector
        else:
            # Wedge doesn't need an inner ring
            v = numpy.vstack([arc.vertices,[(0,0),arc.vertices[0,:],(0,0)]])
            c = numpy.hstack([arc.codes,[connector,connector,Path.CLOSEPOLY]])

        v *= r2
        v += numpy.array(center)
        self._path = Path(v,c)
        self._patch_transform = transforms.IdentityTransform()
    __init__.__doc__ = cbook.dedent(__init__.__doc__) % artist.kwdocd

    def get_path(self):
        return self._path














def readfile(currfile):
    mydatareader=readncnr.datareader()
    mydata=mydatareader.readbuffer(currfile)
    #print mydata.data.keys()
    return mydata

if __name__=='__main__':

    myfilebase='SrFeA0'
    mydirectory=r'C:\srfeas\SrFeAsNi\Ni0p08\May1_2009'
    myend='bt9'
    myfilebaseglob=myfilebase+'*.'+myend
    
    
    
    if 0:
        #Check on lattice constant
        myfilebase='Cspl20'
        mydirectory=r'C:\srfeas\SrFeAsNi\Ni0p08\Ni0p08\2009-03-diffraction'
        myend='bt9'
        myfilebaseglob=myfilebase+'*.'+myend
        myzero=str(0)
        #myzero2=myzero+myzero
        Temp=[]
        clatt=[]
        for i in range(1,29):
            if i<10:
                currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
            else:
                currfile=os.path.join(mydirectory,myfilebase+str(i)+r"."+myend)
            print 'currfile', currfile
            mydata=readfile(currfile)
            T=N.array(mydata.data['temp'])
            a=mydata.metadata['lattice']['a']
            b=mydata.metadata['lattice']['b']
            c=mydata.metadata['lattice']['c']
            Temp.append(T.mean())
            clatt.append(c)
            #mon0=mydata.metadata['count_info']['monitor']
            #mydata=readfile(currfile)
            #print mydata.data.keys()
            #T=N.array(mydata.data['temp'])
            #mon=mydata.metadata['count_info']['monitor']
            #I=N.array(mydata.data['counts'])*mon0/mon
            #Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        pylab.plot(Temp,clatt,'s')
        pylab.show() 
        sys.exit()

    if 1:
        #Summary plots
        myfilebase='Cspl20'
        mydirectory=r'C:\srfeas\SrFeAsNi\Ni0p08\Ni0p08\2009-03-diffraction'
        myend='bt9'
        myfilebaseglob=myfilebase+'*.'+myend
        myzero=str(0)
        #myzero2=myzero+myzero
        Temp=[]
        clatt=[]
        fig=pylab.figure()
        for i in range(1,29):
            if i<10:
                currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
            else:
                currfile=os.path.join(mydirectory,myfilebase+str(i)+r"."+myend)
            print 'currfile', currfile
            mydata=readfile(currfile)
            if i==1:
                mon0=mydata.metadata['count_info']['monitor']
            mydata=readfile(currfile)
            mon=mydata.metadata['count_info']['monitor']
            I=N.array(mydata.data['counts'])*mon0/mon
            Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
            T=N.array(mydata.data['temp'])
            a3=N.array(mydata.data['a3'])
            a=mydata.metadata['lattice']['a']
            b=mydata.metadata['lattice']['b']
            c=mydata.metadata['lattice']['c']
            Temp.append(T.mean())
            clatt.append(c)
            #print mydata.data.keys()
            #T=N.array(mydata.data['temp'])
            #
            #pylab.subplot(6,5,i)
            ax2 = fig.add_subplot(6,5,i)
            pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
            pylab.text(100.6,1000,str(T.mean()),fontsize=12)
            pylab.xlim((100.4,101.8))
            pylab.ylim((0,1400))
            ax2.xaxis.set_major_formatter(pylab.NullFormatter())
            ax2.xaxis.set_major_locator(pylab.NullLocator())
            ax2.yaxis.set_major_formatter(pylab.NullFormatter())
            ax2.yaxis.set_major_locator(pylab.NullLocator())
                #ax2 = fig.add_subplot(2,2,1)
    #pylab.subplots_adjust(hspace=0.6,wspace=0.3)
    #ax2.set_ylim(oymin[center], oymax[center])
    #ax2.yaxis.tick_right()
    #ax2.yaxis.set_label_position('right')
    
            
            
            
        #pylab.plot(Temp,clatt,'s')
        pylab.show()
        sys.exit()

    if 0:
        #Plot the preliminary T- scans
        i=1
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        I=N.array(mydata.data['counts'])
        mon0=mydata.metadata['count_info']['monitor']
        Ierr=N.sqrt(I)
        pylab.subplot(411)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        
    
        
        i=2
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))
        #pylab.annotate('(a)',(140,4000),fonstsize=20)
        #fig1=pylab.gca()
        
        pylab.text(140,3500,'(a)',fontsize=20)
     
        i=4
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        pylab.subplot(412)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))
        pylab.text(140,3500,'(b)',fontsize=20)
     
     
        i=5
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        pylab.subplot(413)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))
     
      
        i=6
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(413)
        cut=29
        pylab.errorbar(T[0:cut],I[0:cut],Ierr[0:cut],marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        pylab.errorbar(T[cut::],I[cut::],Ierr[cut::],marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))
        pylab.text(140,3500,'(c)',fontsize=20)
    
        
        
        
        
    
        i=7
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        pylab.subplot(414)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))
     
        i=8
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(414)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))        
        
        i=9
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T=N.array(mydata.data['temp'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        pylab.errorbar(T,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        pylab.xlabel('T (K)')
        pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((5,160))
        pylab.ylim((2800,4200))    
        pylab.text(140,3500,'(d)',fontsize=20)

    if 0:
        #a3 shifts for the 0 0 4 
        i=13
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T13=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        mon0=mon
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p13=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        #pylab.xlim((103,106))
        
        
        i=18
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T18=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p18=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='green',mec='green',ecolor=None)

        i=19
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T19=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p19=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='black',mec='black',ecolor=None)

        
        

        i=22
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T22=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p22=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        
        
        
        i=23
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T23=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p23=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='yellow',mec='yellow',ecolor=None)

  
        
        
        
    if 0:
        #a3 shifts for the 1 0 3 magnetic peak 
        i=25
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T25=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        mon0=mon
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p25=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        #pylab.xlim((103,106))
        
        
        i=28
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T28=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p28=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='green',mec='green',ecolor=None)

        i=29
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T29=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p29=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='black',mec='black',ecolor=None)

        
        

        i=30
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T30=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p30=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        
        
        
        i=31
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T31=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p31=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='yellow',mec='yellow',ecolor=None)

        
        i=34
        myfilebase='tmags'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T34=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p34=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='orange',mec='orange',ecolor=None)



 
        
    if 1:
        #a3 shifts for the 2 0 6 
        i=34
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T34=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        mon0=mon
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p34=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='blue',mec='blue',ecolor=None)
        #pylab.xlabel('T (K)')
        #pylab.ylabel('Counts //'+str(mon0))
        #pylab.xlim((103,106))
        
        
        i=35
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T35=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p35=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='green',mec='green',ecolor=None)

        i=38
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T38=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p38=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='black',mec='black',ecolor=None)

        
        

        i=39
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T39=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p39=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='red',mec='red',ecolor=None)
        
        
        
        i=40
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T40=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p40=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='yellow',mec='yellow',ecolor=None)

 
        i=41
        myfilebase='tnucs'
        myzero=str(0)
        currfile=os.path.join(mydirectory,myfilebase+myzero+str(i)+r"."+myend)
        print 'currfile', currfile
        mydata=readfile(currfile)
        print mydata.data.keys()
        T41=N.array(mydata.data['temp'])
        a3=N.array(mydata.data['a3'])
        mon=mydata.metadata['count_info']['monitor']
        I=N.array(mydata.data['counts'])*mon0/mon
        Ierr=N.sqrt(N.array(mydata.data['counts']))*mon0/mon
        #pylab.subplot(312)
        p41=pylab.errorbar(a3,I,Ierr,marker='s',linestyle='None',mfc='orange',mec='orange',ecolor=None)

        
        
        
        
        pylab.xlabel('A3 (K)')
        pylab.ylabel('Counts //'+str(mon0))
        pylab.xlim((80,83))
        pylab.ylim((0,9000))
        #print str(T13.mean()), str(T18.mean())
        legstr=[str(T34.mean()), str(T35.mean()), str(T38.mean()),  str(T39.mean()), str(T40.mean())
                ,str(T41.mean())
                
                ]
        legart=[p34[0],p35[0],p38[0],p39[0],p40[0],p41[0]]
        pylab.legend(legart,legstr)
        #pylab.ylim((2800,4200))
        #pylab.annotate('(a)',(140,4000),fonstsize=20)
        #fig1=pylab.gca()
        
        #pylab.text(140,3500,'(a)',fontsize=20)
        
        
        
    pylab.show()
    

    
    #pylab.subplot(311)
    
    #file_range=(43,69)
    
    #mydirectory=r'C:\srfeas\SrFeAsNi\Ni0p08\2009-04-diffraction'
    #file_range=(35,51)
    #myfilebase='SrFeA0'
    #flist=[]
    #for i in range(file_range[0],file_range[1]):
    #    currfile=os.path.join(mydirectory,myfilebase+str(i)+r"."+myend)
    #    #print 'currfile',currfile
    #    flist.append(currfile)