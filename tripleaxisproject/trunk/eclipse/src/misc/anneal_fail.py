import numpy as N
import scipy.optimize.anneal as anneal

def fcn(p):
    x,y=p   
    result=x**2+y**2
    return result


if __name__=="__main__":
    p0=N.array([3,3],'d')
    lowerm=[1,1]
    upperm=[4,4]
    myschedule='fast'
    p0,jmin=anneal(fcn,p0,\
                  schedule=myschedule,lower=lowerm,upper=upperm,\
                  maxeval=None, maxaccept=None,dwell=10,maxiter=600,T0=10000) 
    print 'p0',p0,'jmin',jmin

 