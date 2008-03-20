import numpy as N
import pylab


if __name__=='__main__':
    #1 1 0
    x1=N.array([1])
    y1=N.array([1])
    pylab.subplot(2,1,1)
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
    #-1 -1 0
    x1=N.array([-1])
    y1=N.array([-1])
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
    #-1 2 0
    x1=N.array([-1])
    y1=N.array([2])
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
    #1 -2 0
    x1=N.array([1])
    y1=N.array([-2])
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
    #2 -1 0
    x1=N.array([2])
    y1=N.array([-1])
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
    #-2 1 0
    x1=N.array([-2])
    y1=N.array([1])
    pylab.plot(x1,y1,'bo',markersize=20,markerfacecolor='blue',markeredgecolor='blue')
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



    s=r'$\delta$ 2$\bar{\delta}$ 0'
    pylab.text(1.2,-2,s,fontsize=20)

    s=r'$2\delta$ $\bar{\delta}$ 0'
    pylab.text(1.8,-0.7,s,fontsize=20)

    s=r'$\delta$ $\delta$ 0'
    pylab.text(1,1.2,s,fontsize=20)

    s=r'$\bar{\delta}$ $\bar{\delta}$ 0'
    pylab.text(-1,-1.6,s,fontsize=20)

    s=r'$2\bar{\delta}$ $\delta$ 0'
    pylab.text(-2,1.2,s,fontsize=20)

    s=r'$\bar{\delta}$ 2$\delta$ 0'
    pylab.text(-.7,2.0,s,fontsize=20)

    pylab.axis([-3.5,3.5,-3.5,3.5])
    pylab.show()



