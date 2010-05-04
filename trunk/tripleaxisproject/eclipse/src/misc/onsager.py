import numpy as N
import pylab


if __name__=="__main__":
    jx=1
    jy=100*jx
    T=N.arange(0,10,.02)
    m=(1-1.0/(N.sinh(jx/T)*N.sinh(jy/T))**2)**(.125)
    pylab.plot(T,m,'s')
    pylab.show()