import numpy as N
import lattice_calculator
import rescalc
import prefdemo
#import sqwdemo
#import smademo
import prefbifeo3
import smabifeo3
import convres
import convres_sma

import  pylab
import scipy.sandbox.delaunay as D
import matplotlib.numerix.ma as ma
from matplotlib.ticker import NullFormatter, MultipleLocator
from scipy.signal.signaltools import convolve2d
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator


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


    ax=fig.add_subplot(1,1,nfig)
    pc=ax.pcolormesh(xa,ya,zima,shading='interp',cmap=cmap)  # working good!
#    pc=ax.imshow(zima,interpolation='bilinear',cmap=cmap)
    ##pc.set_clim(0.0,660.0)



    if colorflag:
        #g=pylab.colorbar(pc,ticks=N.arange(0,675,100))
        g=pylab.colorbar(pc)
        print g
        #g.ticks=None
        #gax.yaxis.set_major_locator(MultipleLocator(40))
        #g.ticks(N.array([0,20,40,60,80]))

    return ax,g


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







if __name__=="__main__":
    if 1:
        a=N.array([3.97296],'d')
        b=N.array([2.79045],'d')
        c=N.array([13.86],'d')
        alpha=N.array([90],'d')
        beta=N.array([90],'d')
        gamma=N.array([90],'d')
 #       orient1=N.array([[0,1,1]],'d')
        orient1=N.array([[1,0,0]],'d')
        orient2=N.array([[0,-1,0]],'d')
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)
        H=N.array([.5,1],'d');K=N.array([.5,1.2],'d');L=N.array([0],'d');W=N.array([0],'d')
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
        xmin=.45
        xmax=.55
        ymin=.45
        ymax=.55
        #.485, .515
        ht=N.linspace(xmin,xmax,5)
        kt=N.linspace(ymin,ymax,5)
        h,k=N.meshgrid(ht,kt)
        H=h.flatten()
        K=k.flatten()
        #print 'H shape ',H.shape
        #print 'K shape ',K.shape
        #Parameter values for the cross section
        p=N.zeros((4,),'float64')
        #p=[hcenter,kcenter,lcenter,corr,Intensity,background]
        p=[0.5,0.5,0,2.0e1,1000,0]
        myrescal=rescalc.rescalculator(mylattice)
        newinput=lattice_calculator.CleanArgs(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orient1=orient1,orient2=orient2,\
                            H=H,K=K,L=L,W=W,setup=setup)
        mylattice=lattice_calculator.lattice(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                        beta=newinput['beta'],gamma=newinput['gamma'],orient1=newinput['orient1'],\
                        orient2=newinput['orient2'])
        myrescal.__init__(mylattice)
        H=newinput['H']
        K=newinput['K']
        L=newinput['L']
        W=newinput['W']
        setup=newinput['setup']
        #R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
        R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
        #print 'RMS'
        #print RMS.transpose()[0]

        #myrescal.ResPlot(H, K, L, W, setup)

        (prefactor,background)=prefbifeo3.PrefDemo(H,K,L,W,myrescal,p)
        #print 'prefactor ',prefactor
        #print 'background ',background
        #sqw=sqwdemo.SqwDemo(H,K,L,W,p)
        disp,myint,WL=smabifeo3.SMADemo(H,K,L,p)
        #print disp
        #print myint
        #print WL
        #print myint.shape
        ac=[5,0]
        conv_sma=convres_sma.ConvResSMA(smabifeo3.SMADemo,prefbifeo3.PrefDemo,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=ac)
        print 'conv '
        print conv_sma.shape
        print conv_sma.max()
        print myint.max()
        exit()
        #print 'sqw ', sqw
        xd,yd,zd=prep_data2(H,K,conv_sma[0,:]);
        if 1:

            fig=pylab.figure(figsize=(8,8))
            ylim=(ymin,ymax)
            xlim=(xmin,xmax)
            #ylabel='E (meV)'
            #xlabel=r'Q$ \ \ (\AA^{-1}$)'
            ylabel='(1 1 0)'
            xlabel='(1 -1 -2)'


        if 1:
            ax,g=plot_data(xd,yd,zd,fig,1,colorflag=True)
            #ax.text(.98,.20,'E=0 KV',fontsize=14,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes,color='white')
            ax.set_ylabel(ylabel)
            ax.set_xlabel(xlabel)
            ax.xaxis.set_major_formatter(NullFormatter())
            ax.set_ylim(ylim); ax.set_xlim(xlim)
            #g.ax.ticks=N.arange(0,100,20)

        if 1:
            print 'gca ', fig.gca()
            for im in fig.gca().get_images():
                print im
                ###im.set_clim(0.0,660.0) #should reset
            pylab.show()



        exit()
        ac=[5,0]
        #conv=convres.ConvRes(sqwdemo.SqwDemo,prefdemo.PrefDemo,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=ac)
        #print conv
        conv_sma=convres_sma.ConvResSMA(smabifeo3.SMADemo,prefbifeo3.PrefDemo,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=ac)
        print 'conv'
        #print conv_sma
        print 'shape ',conv_sma.shape
        #print 'H shape ', H.shape
        #print 'W shape ', W.shape
        #R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
        #myrescal.ResPlot(H, K, L, W, setup)
        #print 'RMS'
        #print RMS.transpose()[0]
        #print myrescal.calc_correction(H,K,L,W,setup,qscan=[[1,1,0],[1,1,0]])
        #print myrescal.CalcWidths(H,K,L,W,setup)