import numpy as N
import pylab
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.ticker import NullFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator,NullLocator


if __name__=='__main__':
    #layer 1 z=0.125
    x=[0.5,0.5,0.5,0.5,1.5,1.5,1.5,1.5]
    y=[0.25,0.75,1.25,1.75,0.25,0.75,1.25,1.75]
    theta=0
    a=0.1
    fig=pylab.figure()
    ax = fig.add_subplot(111)
    x0y0=N.array([0.0,0.0])
    r=Rectangle(x0y0,height=1,width=1,fill=True,facecolor='white',edgecolor='black',linewidth=3.0)
    ax.add_artist(r)
    for i in range(len(x)):
        x0y0=N.array([x[i],y[i]])
        e=Ellipse(x0y0,width=2*a,height=2*a,angle=theta)
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        #e.set_alpha(0.5)
        e.set_facecolor('red')
    ax.set_xlim(-.1,2.1)
    ax.set_ylim(-.1,2.1)
    ax.yaxis.set_major_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(NullLocator())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
   
    
    pylab.box(on=None)
    pylab.show()