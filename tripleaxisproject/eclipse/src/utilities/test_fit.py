from __future__ import division
import copy
import numpy as N
pi=N.pi
from mpfit import mpfit
import matplotlib
import pylab

def gauss(p,x):
    #Area center width Bak

    #p=[p0,p1,p2,p3]


    x0=p[1]
    width=p[2]
    sigma=width/2/N.sqrt(2*N.log(2))
    area=N.abs(p[0])/N.sqrt(2*pi)/sigma
    background=N.abs(p[3])
    y=background+area*N.exp(-(0.5*(x-x0)*(x-x0)/sigma/sigma))
    return y

def chisq(p,a3,I,Ierr):
    Icalc=gauss(p,a3)
    #print I.shape
    #print Ierr.shape
    #print a3.shape
    #print Icalc.shape
    Ierr_temp=copy.deepcopy(Ierr)
    zero_loc=N.where(Ierr==0)[0]
    if len(zero_loc)!=0:
        Ierr_temp[zero_loc]=1.0
    chi=((I-Icalc)/Ierr_temp)**2    
    return chi.sum()/(len(I)-len(p))


def myfunct_res(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = gauss(p, x)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    Ierr_temp=copy.deepcopy(err)
    zero_loc=N.where(err==0)[0]
    if len(zero_loc)!=0:
        Ierr_temp[zero_loc]=1.0
    return [status, (y-model)/Ierr_temp]


def prepfit(x,y,yerr,area=3.0,center=2.0,width=2.0,Bak=5.0):
    p0=N.array([area,center,width,Bak],dtype='float64')  #initial conditions
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)): 
        parinfo[i]['value']=p0[i]
    fa = {'x':x, 'y':y, 'err':yerr}
    parinfo[1]['fixed']=0
    parinfo[2]['fixed']=0
    m = mpfit(myfunct_res, p0, parinfo=parinfo,functkw=fa)
    if (m.status <= 0): 
        print 'error message = ', m.errmsg
    params=m.params
    pfit=params
    perror=m.perror
    print 'pfit',pfit
    #chisqr=(myfunct_res(m.params, x=th, y=counts, err=counts_err)[1]**2).sum()
    chisqr=chisq(pfit,x,y,yerr)
    dof=m.dof
    Icalc=gauss(pfit,x)
    return Icalc


if __name__=="__main__":
    x=N.arange(-5,5,.05)
    area=150.0
    center=0.0
    width=2.0
    background=4.0
    p=N.array([area, center, width, background])
    y=gauss(p,x)
    yerr=N.sqrt(y)
    Icalc=prepfit(x,y,yerr,area=area,center=center, width=width, Bak=background)
    if 1:
        pylab.errorbar(x,y,yerr,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
        pylab.plot(x,Icalc)
        pylab.show()
    
    