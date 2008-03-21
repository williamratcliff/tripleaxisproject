import numpy as N
import pylab
import math
import unittest
from matplotlib.patches import Ellipse
import lattice_calculator
import readncnr2 as readncnr
import scriptutil as SU
import re
import simple_combine
import sys
sys.path.append(r'c:\tripleaxisproject2\polarization')
import polcorrect_lib
from matplotlib.ticker import NullFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
import scipy.sandbox.delaunay as D
#import numpy.core.ma as ma
import matplotlib.numerix.ma as ma
from scipy.signal.signaltools import convolve2d
import readicp



eps=1e-3
pi=N.pi


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

def readfiles(mydirectory,myfilebase,myend):
    myfilebaseglob=myfilebase+'*.'+myend
    print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    Qx=N.array([])
    Qy=N.array([])
    Qz=N.array([])
    Counts=N.array([])
    for currfile in flist:
        print currfile
        mydata=mydatareader.readbuffer(currfile)
        Qx=N.concatenate((Qx,N.array(mydata.data['qx'])))
        Qy=N.concatenate((Qy,N.array(mydata.data['qy'])))
        Qz=N.concatenate((Qz,N.array(mydata.data['qz'])))
        Counts=N.concatenate((Counts,N.array(mydata.data['counts'])))


def plot_nodes(tri):
    for nodes in tri.triangle_nodes:
        D.fill(x[nodes],y[nodes],'b')
    pylab.show()

def plot_data(xa,ya,za,fig,nfig,colorflag=False):

    cmap = pylab.cm.jet
    cmap.set_bad('w', 1.0)
    myfilter=N.array([[0.1,0.2,0.1],[0.2,0.8,0.2],[0.1,0.2,0.1]],'d') /2.0
    zout=convolve2d(za,myfilter,mode='same')
    zima = ma.masked_where(N.isnan(zout),zout)


    ax=fig.add_subplot(3,2,nfig)
    pc=ax.pcolormesh(xa,ya,zima,shading='interp',cmap=cmap)  # working good!
#    pc=ax.imshow(zima,interpolation='bilinear',cmap=cmap)
    pc.set_clim(0.0,660.0)



    if colorflag:
        g=pylab.colorbar(pc,ticks=N.arange(0,675,100))
        print g
        #g.ticks=None
        #gax.yaxis.set_major_locator(MultipleLocator(40))
        #g.ticks(N.array([0,20,40,60,80]))

    return ax,g

def prep_data(filename):
#    Data=pylab.load(r'c:\resolution_stuff\1p4K.iexy')
    Data=pylab.load(filename)
    xt=Data[:,2]
    yt=Data[:,3]
    zorigt=Data[:,0]
    x=xt[:,zorigt>0.0]
    y=yt[:,zorigt>0.0]
    z=zorigt[:,zorigt>0.0]
#    zorig=ma.array(zorigt)
    print 'reached'
    threshold=0.0;
#    print zorigt < threshold
#    print N.isnan(zorigt)
#    z = ma.masked_where(zorigt < threshold , zorigt)
    print 'where masked ', z.shape
#should be commented out--just for testing
##    x = pylab.randn(Nu)/aspect
##    y = pylab.randn(Nu)
##    z = pylab.rand(Nu)
##    print x.shape
##    print y.shape
    # Grid
    xi, yi = N.mgrid[-5:5:100j,-5:5:100j]
    xi,yi=N.mgrid[x.min():x.max():.05,y.min():y.max():.05]
    # triangulate data
    tri = D.Triangulation(x,y)
    print 'before interpolator'
    # interpolate data
    interp = tri.nn_interpolator(z)
    print 'interpolator reached'
    zi = interp(xi,yi)
    # or, all in one line
    #    zi = Triangulation(x,y).nn_interpolator(z)(xi,yi)
#    return x,y,z
    return xi,yi,zi





def prep_data2(xt,yt,zorigt):
#    Data=pylab.load(r'c:\resolution_stuff\1p4K.iexy')
    #Data=pylab.load(filename)
    #xt=Data[:,2]
    #yt=Data[:,3]
    #zorigt=Data[:,0]
    x=xt[:,zorigt>0.0]
    y=yt[:,zorigt>0.0]
    z=zorigt[:,zorigt>0.0]
#    zorig=ma.array(zorigt)
    print 'reached'
    threshold=0.0;
#    print zorigt < threshold
#    print N.isnan(zorigt)
#    z = ma.masked_where(zorigt < threshold , zorigt)
    print 'where masked ', z.shape
#should be commented out--just for testing
##    x = pylab.randn(Nu)/aspect
##    y = pylab.randn(Nu)
##    z = pylab.rand(Nu)
##    print x.shape
##    print y.shape
    # Grid
    xi, yi = N.mgrid[-5:5:100j,-5:5:100j]
    xi,yi=N.mgrid[x.min():x.max():.001,y.min():y.max():.001]
    # triangulate data
    tri = D.Triangulation(x,y)
    print 'before interpolator'
    # interpolate data
    interp = tri.nn_interpolator(z)
    print 'interpolator reached'
    zi = interp(xi,yi)
    # or, all in one line
    #    zi = Triangulation(x,y).nn_interpolator(z)(xi,yi)
#    return x,y,z
    return xi,yi,zi



def readmeshfiles(mydirectory,myfilebase,myend):
    myfilebaseglob=myfilebase+'*.'+myend
    print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    mydatareader=readicp.datareader()
    Qx=N.array([])
    Qy=N.array([])
    Qz=N.array([])
    Counts=N.array([])
    for currfile in flist:
        print currfile
        mydata=mydatareader.readbuffer(currfile)
        Qx=N.concatenate((Qx,N.array(mydata.data['Qx'])))
        Qy=N.concatenate((Qy,N.array(mydata.data['Qy'])))
        Qz=N.concatenate((Qz,N.array(mydata.data['Qz'])))
        Counts=N.concatenate((Counts,N.array(mydata.data['Counts'])))
    xa,ya,za=prep_data2(Qx,Qy,Counts);
    return xa,ya,za

def strangezone(fig):
    Nu = 10000
    aspect = 1.0
    mydirectory=r'c:\bifeo3xtal\dec7_2007'
    myfilebase='cmesh'
    myend='bt9'
    xc,yc,zc=readmeshfiles(mydirectory,'cmesh',myend) #Rm temp
    xd,yd,zd=readmeshfiles(mydirectory,'dmesh',myend) #0
    xe,ye,ze=readmeshfiles(mydirectory,'emesh',myend) #-1.3
    xf,yf,zf=readmeshfiles(mydirectory,'fmesh',myend) #0
    xg,yg,zg=readmeshfiles(mydirectory,'gmesh',myend) #1.3
    xh,yh,zh=readmeshfiles(mydirectory,'hmesh',myend) #0
    xi,yi,zi=readmeshfiles(mydirectory,'imesh',myend) #-1.3
    xj,yj,zj=readmeshfiles(mydirectory,'jmesh',myend) #0
    print 'matplotlib'

    if 1:

        #fig=pylab.figure(figsize=(8,8))
        ylim=(.485,.515)
        xlim=(.485,.515)
        #ylabel='E (meV)'
        #xlabel=r'Q$ \ \ (\AA^{-1}$)'
        ylabel='(1 1 0)'
        xlabel='(1 -1 -2)'


    if 1:
        ax,g=plot_data(xd,yd,zd,fig,4,colorflag=True)
        #ax.text(.98,.20,'E=0 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
        #g.ax.ticks=N.arange(0,100,20)

    if 0:
        ax,g=plot_data(xe,ye,ze,fig,2,colorflag=True)
        ax.text(.98,.20,'E=-1.3 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)

    if 0:
        ax,g=plot_data(xf,yf,zf,fig,3,colorflag=True)
        ax.text(.98,.20,'0 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
    if 0:
        ax,g=plot_data(xg,yg,zg,fig,4,colorflag=True)
        ax.text(.98,.20,'1.3 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
    if 0:
        ax,g=plot_data(xh,yh,zh,fig,5,colorflag=True)
        ax.text(.98,.20,'0 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
        fmt = FormatStrFormatter('%0.3g')  # or whatever
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(MaxNLocator(5))
    if 0:
        ax,g=plot_data(xi,yi,zi,fig,6,colorflag=True)
        ax.text(.98,.20,'-1.3 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.set_xlabel(xlabel)
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
        fmt = FormatStrFormatter('%0.3g')  # or whatever
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(MaxNLocator(5))

    if 0:
        ax,g=plot_data(xf,yf,zf,fig,6,colorflag=True)
        ax.text(.98,.20,'0 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
        ax.set_xlabel(xlabel)
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.set_ylim(ylim); ax.set_xlim(xlim)
        #fmt = FormatStrFormatter('%1.4g')  # or whatever
        #ax.yaxis.set_major_formatter(fmt)


    if 1:
        print 'gca ', fig.gca()
        for im in fig.gca().get_images():
            print im
            im.set_clim(0.0,660.0)
        #pylab.show()
    if 0:
        print 'saving'
        pylab.savefig(r'c:\sqltest\demo.pdf',dpi=150)
        print 'saved'
    if 0:
        pylab.show()



def cartoon():
    if 1:
        #h k 0
        ax=pylab.subplot(3,4,1)
        #1 1 0
        x1=N.array([1])
        y1=N.array([1])

        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        #-1 -1 0
        x1=N.array([-1])
        y1=N.array([-1])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        #-1 2 0
        x1=N.array([-1])
        y1=N.array([2])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        #1 -2 0
        x1=N.array([1])
        y1=N.array([-2])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        #2 -1 0
        x1=N.array([2])
        y1=N.array([-1])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        #-2 1 0
        x1=N.array([-2])
        y1=N.array([1])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        pylab.xlabel('[1 0 0]')
        pylab.ylabel('[0 1 0]')
        #draw vertical line
        #x1=[0,0]
        #y1=[-2.5,2.5]
        #pylab.plot(x1,y1,linewidth=3.0)
        #draw horizontal line
        #x1=[-2.5,2.5]
        #y1=[0,0]
        #pylab.plot(x1,y1,linewidth=3.0,color='blue')

        #draw line 1 1 0
        x1=[1,-1]
        y1=[1,-1]
        pylab.plot(x1,y1,linewidth=3.0,color='blue')

        #draw line 1 2 0
        x1=[1,-1]
        y1=[-2,2]
        pylab.plot(x1,y1,linewidth=3.0,color='blue')

        #draw line 1 1 0
        x1=[-2,2]
        y1=[1,-1]
        pylab.plot(x1,y1,linewidth=3.0,color='blue')



        s=r'$\delta 2\bar{\delta} 0$'
        pylab.text(1.2,-2,s,fontsize=10)

        s=r'$2\delta \bar{\delta} 0$'
        pylab.text(1.8,-0.7,s,fontsize=10)

        s=r'$\delta \delta 0$'
        pylab.text(1,1.2,s,fontsize=10)

        s=r'$\bar{\delta} \bar{\delta} 0$'
        pylab.text(-1,-1.6,s,fontsize=10)

        s=r'$2\bar{\delta} \delta 0$'
        pylab.text(-2,1.2,s,fontsize=10)

        s=r'$\bar{\delta} 2\delta 0$'
        pylab.text(-.7,2.0,s,fontsize=10)

        pylab.axis([-3.5,3.5,-3.5,3.5])
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())

    if 1:
        #h h l
        ax=pylab.subplot(3,4,2)

        # 1 1 0
        x1=N.array([4])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        s=r'$\delta \delta 0$'
        pylab.text(4.8,0.0,s,fontsize=10)

        # -1 -1 0
        x1=N.array([-4])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='blue',markeredgecolor='blue')
        s=r'$\bar{\delta} \bar{\delta} 0$'
        pylab.text(-6.8,0.0,s,fontsize=10)

        # -1 2 0
        x1=N.array([2])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=5,markerfacecolor='white')
        s=r'$\bar{\delta} 2\delta 0$'
        pylab.text(0.5,0.4,s,fontsize=10)

        # 2 -1 0
        #x1=N.array([2])
        #y1=N.array([0])
        #pylab.plot(x1,y1,'bo',markersize=10)
        s=r'$2\delta \bar{\delta} 0$'
        pylab.text(0.5,-1.4,s,fontsize=10)


        # -2 1 0
        x1=N.array([-2])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=5,markerfacecolor='white')
        s=r'$2\bar{\delta} \delta 0$'
        pylab.text(-3.5,0.4,s,fontsize=10)

        # 1 -2 0
        #x1=N.array([2])
        #y1=N.array([0])
        #pylab.plot(x1,y1,'bo',markersize=10)
        s=r'$\delta 2\bar{\delta} 0$'
        pylab.text(-3.5,-1.4,s,fontsize=10)



        pylab.axis([-8.0,8.0,-3.5,3.5])
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())

        pylab.xlabel('[1 1 0]')
        pylab.ylabel('[0 0 1]')

        #ax=pylab.subplot(4,2,2)
    if 1:
        #h k h-k
        ax=pylab.subplot(3,2,2)
        # 1 1 0
        x1=N.array([4])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
        s=r'$\delta \delta 0$'
        pylab.text(4.8,0.0,s,fontsize=20)

        # -1 -1 0
        x1=N.array([-4])
        y1=N.array([0])
        pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
        s=r'$\bar{\delta} \bar{\delta} 0$'
        pylab.text(-6.8,0.0,s,fontsize=20)

        # -1 2 0
        x1=N.array([2])
        y1=N.array([4])
        pylab.plot(x1,y1,'bo',markersize=5,markerfacecolor='white')
        s=r'$\bar{\delta} 2\delta 0$'
        pylab.text(0.5,4.4,s,fontsize=20)


        # -2 1 0
        x1=N.array([-2])
        y1=N.array([4])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='white')
        s=r'$2\bar{\delta} \delta 0$'
        pylab.text(-3.0,4.4,s,fontsize=20)


        # 2 -1 0
        x1=N.array([2])
        y1=N.array([-4])
        pylab.plot(x1,y1,'bo',markersize=10,markerfacecolor='gray')
        s=r'$2\delta \bar{\delta} 0$'
        pylab.text(0.5,-6.8,s,fontsize=20)


        # 1 -2 0
        x1=N.array([-2])
        y1=N.array([-4])
        pylab.plot(x1,y1,'bo',markersize=5,markerfacecolor='gray')
        s=r'$\delta 2\bar{\delta} 0$'
        pylab.text(-3.0,-6.8,s,fontsize=20)


        pylab.axis([-10.5,10.5,-10.5,10.5])
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())

        pylab.xlabel('[1 1 0]')
        pylab.ylabel('[-1 1 -2]')

    if 0:
        pylab.subplot(3,2,3)
        pylab.subplot(3,2,4)
        pylab.subplot(3,2,5)
        pylab.subplot(3,2,6)
        pylab.subplots_adjust(wspace=0.5,hspace=0.6)
    return


if __name__=="__main__":
    if 1:
        fig=pylab.figure(figsize=(8,8))
        a=N.array([5.5436],'d')
        b=N.array([5.5436],'d')
        c=N.array([13.8844],'d')
        alpha=N.array([90],'d')
        beta=N.array([90],'d')
        gamma=N.array([120],'d')
 #       orient1=N.array([[0,1,1]],'d')
        orient1=N.array([[0,0,1]],'d')
        orient2=N.array([[1,1,0]],'d')
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)

        EXP={}
        EXP['ana']={}
        EXP['ana']['tau']='pg(002)'
        EXP['mono']={}
        EXP['mono']['tau']='pg(002)';
        EXP['ana']['mosaic']=25
        EXP['mono']['mosaic']=25
        EXP['sample']={}
        EXP['sample']['mosaic']=25
        EXP['sample']['vmosaic']=25
        EXP['hcol']=N.array([40, 40, 40, 80],'d')
        EXP['vcol']=N.array([120, 120, 120, 120],'d')
        EXP['infix']=-1 #positive for fixed incident energy
        EXP['efixed']=14.7
        EXP['method']=0
        setup=[EXP]

        mydirectory=r'c:\BiFeO3xtal\Oct9_2007'
        myfilebase='bfo_spinflip_he3flip51585'
        myfilebase2='bfo_spinflip51579'
        myfilebase3='bfo_spinflip51583'
        myend='bt7'
        data={}
        #pm
        myfilebaseglob=myfilebase+'*.'+myend
        print myfilebaseglob
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilestr=flist[0]
        Counts=N.array([],'float64')
        monlist=[]
        timestamp={}
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        S1=N.array(mydata.data['a3'][1:])
        S2=N.array(mydata.data['a4'][1:])
        A2=N.array(mydata.data['a5'][1:])*2
        M2=N.array(mydata.data['a2'][1:])
        Counts_pm=N.concatenate((Counts,N.array(mydata.data['detector'][1:])))
        H_pm,K_pm,L_pm,E_pm,Q_pm,Ei_pm,Ef_pm=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
        monlist.append(mydata.data['monitor'][0])
        Counts_pm=Counts_pm[:-1]
        H_pm=H_pm[:-1]
        K_pm=K_pm[:-1]
        L_pm=L_pm[:-1]
        E_pm=E_pm[:-1]
        Q_pm=Q_pm[:-1]
        Ei_pm=Ei_pm[:-1]
        Ef_pm=Ef_pm[:-1]
        timestamp['pm']=N.array(mydata.data['timestamp'][1:])[:-1]
        #print H_pm
        #begin mp
        myfilebaseglob=myfilebase2+'*.'+myend
        print myfilebaseglob
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilebaseglob=myfilebase3+'*.'+myend
        flist2 = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))

        myfilestr=flist[0]
        Counts=N.array([],'float64')
        mydatareader=readncnr.datareader()
        mydata1=mydatareader.readbuffer(myfilestr)
        S1=N.array(mydata1.data['a3'])[1:]
        S2=N.array(mydata1.data['a4'])[1:]
        A2=N.array(mydata1.data['a5'])[1:]*2
        M2=N.array(mydata1.data['a2'])[1:]
        Counts_mp=N.array(mydata1.data['detector'])[1:]
        timestamp_mp=N.array(mydata1.data['timestamp'])[1:]
        #H_mp,K_mp,L_mp,E_mp,Q_mp,Ei_mp,Ef_mp=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
        monlist.append(mydata1.data['monitor'][0])
        #print M2.shape
        #print Counts_mp.shape
        #remove potential bad datapoint that bt7 generated in dd mode


        myfilestr=flist2[0]
        mydatareader=readncnr.datareader()
        mydata2=mydatareader.readbuffer(myfilestr)
        S1=N.concatenate((S1,N.array(mydata2.data['a3'])[1:]))
        S2=N.concatenate((S2,N.array(mydata2.data['a4'])[1:]))
        A2=N.concatenate((A2,N.array(mydata2.data['a5'])[1:]*2))
        M2=N.concatenate((M2,N.array(mydata2.data['a2'])[1:]))
        timestamp_mp=N.concatenate((timestamp_mp,N.array(mydata2.data['timestamp'])[1:]))
        timestamp['mp']=timestamp_mp
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)

        Counts_mp=N.concatenate((Counts_mp,N.array(mydata2.data['detector'])[1:]))
        #print 'where'
        H_mp,K_mp,L_mp,E_mp,Q_mp,Ei_mp,Ef_mp=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)


        ylist=[Counts_pm,Counts_mp]
        yerrlist=[N.sqrt(Counts_pm),N.sqrt(Counts_mp)]
        I,Ierr=simple_combine.monitor_normalize(ylist,yerrlist,monlist)
        counts={}
        errs={}
        counts['pm']=I[0]
        counts['mp']=I[1]
        errs['pm']=Ierr[0]
        errs['mp']=Ierr[1]
        cell=mydirectory+'\WilliamOct2007horizCells.txt'
        print cell
        pbflags=polcorrect_lib.PBflags()
        pbflags.MonitorCorrect=0
        pbflags.PolMonitorCorrect=1
        pbflags.MonoSelect=1
        pbflags.Debug=0
        pbflags.SimFlux=0
        pbflags.SimDeviate=0
        pbflags.NoNegativeCS=0
        pbflags.HalfPolarized=0
        pbflags.CountsEnable=[0,0,1,1]
        pbflags.CountsAdd1=[0,0,0,0]
        pbflags.CountsAdd2=[0,0,0,0]
        pbflags.Sconstrain=[1,1,0,0]
        pbflags.Spp=[0,0,0,0]
        pbflags.Smm=[0,0,0,0]
        pbflags.Spm=[0,0,0,0]
        pbflags.Smp=[0,0,0,0]

        mypolcor=polcorrect_lib.polarization_correct(counts,errs,timestamp,cell,Ei_mp,Ef_mp)
        corrected_counts=mypolcor.correct(pbflags)

        #print timestamp['mp'].shape
        #print timestamp['pm'].shape
        #print S1.shape
        #print S2.shape
        #print A2.shape
        #print M2.shape
        #print Counts_mp.shape
        #print H_mp.shape
        #print Counts_mp.shape
        #print I[0].shape
        #print I[1].shape



# Next zone
        myend='out'
        mydirectory=r'c:\bifeo3xtal\jan8_2008\9175\data'
        myfilebase='fieldscansplusminusreset53630'
        myfilebase2='fieldscanminusplusreset53631'

        myfilebaseglob=myfilebase+'*.'+myend
        flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
        myfilebaseglob=myfilebase2+'*.'+myend
        flist2 = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))

        myfilestr=flist[0]
        mydatareader=readncnr.datareader()
        mydata_pm=mydatareader.readbuffer(myfilestr)
        q_pm=N.array(mydata_pm.data['qx'])
        counts_pm=N.array(mydata_pm.data['detector_corrected'])
        errs_pm=N.array(mydata_pm.data['detector_errs_corrected'])


        myfilestr=flist2[0]
        mydatareader=readncnr.datareader()
        mydata_mp=mydatareader.readbuffer(myfilestr)
        q_mp=N.array(mydata_mp.data['qx'])
        counts_mp=N.array(mydata_mp.data['detector_corrected'])
        errs_mp=N.array(mydata_mp.data['detector_errs_corrected'])
        cartoon()

        if 1:
            pylab.subplot(3,2,6)
            pylab.errorbar(q_mp,counts_mp,fmt='bo',yerr=errs_mp,linestyle='None')
            pylab.errorbar(q_pm,counts_pm,fmt='ro',yerr=errs_pm,linestyle='None')


        if 1:
            pylab.subplot(3,2,5)
            pylab.errorbar(H_mp,corrected_counts['Smp'],fmt='bo',yerr=corrected_counts['Emp'],linestyle='None')
            pylab.errorbar(H_pm,corrected_counts['Spm'],fmt='ro',yerr=corrected_counts['Epm'],linestyle='None')

        if 0:
            pylab.errorbar(H_mp,data['mp'],fmt='bo',yerr=errs['mp'],linestyle='None')
            pylab.errorbar(H_pm,data['pm'],fmt='ro',yerr=errs['pm'],linestyle='None')
            pylab.show()

        strangezone(fig)
        pylab.subplots_adjust(wspace=0.5,hspace=0.6)
        pylab.show()
        exit()
        if 0:
            for i in range(S1.size):
                s='%f %f %f %f'%(S1[i],S2[i],A2[i],M2[i])
                print s
            exit()
        if 0:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians([M2[0]]),myradians([-163.193]),myradians([29.833]),myradians([A2[0]]),setup)
            s='%f %f %f %f %f %f %f'%(S1[0],S2[0],H[0],K[0],L[0],E[0],Q[0])
            print s
        if 0:
            H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(myradians(M2),myradians(S1),myradians(S2),myradians(A2),setup)
            for i in range(1,S1.size):
                s='%f %f %f %f %f %f %f'%(S1[i],S2[i],H[i],K[i],L[i],E[i],Q[i])
                print s

        if 0:
            pylab.plot(H,L,'bo')
            pylab.xlabel('[1,1,0]')
            pylab.ylabel('[0,0,1]')
            #pylab.axis([-1,1,-4,4])
            pylab.show()

