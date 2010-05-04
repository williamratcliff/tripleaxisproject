import numpy as N
import pylab
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.ticker import NullFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator,NullLocator


def makeboxes(ax):
    x0y0=N.array([0.0,0.0])
    r=Rectangle(x0y0,height=1,width=1,fill=True,facecolor='white',edgecolor='black',linewidth=3.0)
    ax.add_artist(r)
    x0y0=N.array([1.0,0.0])
    r=Rectangle(x0y0,height=1,width=1,fill=True,facecolor='white',edgecolor='black',linewidth=3.0)
    ax.add_artist(r)
    x0y0=N.array([0.0,1.0])
    r=Rectangle(x0y0,height=1,width=1,fill=True,facecolor='white',edgecolor='black',linewidth=3.0)
    ax.add_artist(r)
    x0y0=N.array([1.0,1.0])
    r=Rectangle(x0y0,height=1,width=1,fill=True,facecolor='white',edgecolor='black',linewidth=3.0)
    ax.add_artist(r)

def addellipses(x,y,ax,a=0.1,color='red'):
    for i in range(len(x)):
        x0y0=N.array([x[i],y[i]])
        e=Ellipse(x0y0,width=2*a,height=2*a,angle=theta)
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        #e.set_alpha(0.5)
        e.set_facecolor(color)

def setaxes(ax):
    ax.set_xlim(-.3,2.3)
    ax.set_ylim(-.3,2.3)
    ax.yaxis.set_major_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(NullLocator())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())    
    pylab.box(on=None)

def z0p25(fig):
    ax = fig.add_subplot(221)
    makeboxes(ax)
    #mn
    x=[0.5,0.5,0.5,0.5,
       1.5,1.5,1.5,1.5]
    y=[0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75]   
    addellipses(x,y,ax,a=0.1,color='red')
    #oxygen 
    x=[0.25,0.25,0.25,0.25,0.25,
       0.75,0.75,0.75,0.75,0.75,
       1.25,1.25,1.25,1.25,1.25,
       1.75,1.75,1.75,1.75,1.75]
    y=[0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0]
    addellipses(x,y,ax,a=0.15,color='white')
    setaxes(ax)
    pylab.title('z=0.125')
    #add atom labels
    #pylab.text(0.,0.5,'O')
    #pylab.text(0.5,.25,'Mn')


def z0p375(fig):
    ax = fig.add_subplot(222)
    makeboxes(ax)
    #mn
    x=[0.25,0.25,0.25,
       0.75,0.75,0.75,
       1.25,1.25,1.25,
       1.75,1.75,1.75,
       ]
    y=[0,1,2,
       0,1,2,
       0,1,2,
       0,1,2
       ]   
    addellipses(x,y,ax,a=0.1,color='red')
    #oxygen 
    x=[0,0,0,0,
       0.5,0.5,0.5,0.5,
       1.,1.,1.,1.,
       1.5,1.5,1.5,1.5,
       2.,2.,2.,2.,
       ]
    y=[0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75]
    addellipses(x,y,ax,a=0.15,color='white')
    setaxes(ax)
    pylab.title('z=0.375')

    

def z0p625(fig):
    ax = fig.add_subplot(223)
    makeboxes(ax)
    #mn
    x=[0,0,0,0,
       1,1,1,1,
       2,2,2,2] 
    y=[0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       ]   
    addellipses(x,y,ax,a=0.1,color='red') 
    #oxtgen
    x=[0.25,0.25,0.25,0.25,0.25,
       0.75,0.75,0.75,0.75,0.75,
       1.25,1.25,1.25,1.25,1.25,
       1.75,1.75,1.75,1.75,1.75]
    y=[0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0,
       0,0.5,1,1.5,2.0]
    addellipses(x,y,ax,a=0.15,color='white')
    setaxes(ax)
    pylab.title('z=0.625')



def z0p875(fig):
    ax = fig.add_subplot(224)
    makeboxes(ax)
    #mn
    x=[0.25,0.25,
       0.75,0.75,
       1.25,1.25,
       1.75,1.75,
       ]
    y=[0.5,1.5,
       0.5,1.5,
       0.5,1.5,
       0.5,1.5
       ]   
    addellipses(x,y,ax,a=0.1,color='red')
    #oxygen 
    x=[0,0,0,0,
       0.5,0.5,0.5,0.5,
       1.,1.,1.,1.,
       1.5,1.5,1.5,1.5,
       2.,2.,2.,2.,
       ]
    y=[0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75,
       0.25,0.75,1.25,1.75]
    addellipses(x,y,ax,a=0.15,color='white')
    setaxes(ax)
    pylab.title('z=0.875')





    
if __name__=='__main__':
    #layer 1 z=0.125
    theta=0
    
    fig=pylab.figure()
    z0p25(fig)
    z0p375(fig)
    z0p625(fig)
    z0p875(fig)
    pylab.show()