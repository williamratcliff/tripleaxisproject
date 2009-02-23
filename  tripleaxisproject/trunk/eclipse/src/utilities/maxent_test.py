from __future__ import division
import numpy as N
import dct
import pylab

pi=N.pi

def plotdensity(h,k,l,fq,xstep=0.01,zstep=0.01):
    x=N.arange(0.0,1.0,xstep)
    z=N.arange(0.0,1.0,zstep)
    xn=len(x)
    zn=len(z)
    X,Z=N.meshgrid(x,z)
    fsum=0.0
    P=N.zeros(X.shape)
    delta=.065;
    #%delta=.1;
    delta=.035;
    #outfile='c:\structfactors_density.dat'
    #fid=fopen([outfile],'w');
    for xia in range(xn):
        for zia in range(zn):
            fsum=0
            xi=x[xia]
            zi=z[zia]
            Aj=fq*N.sinc(2*delta*h)*N.sinc(2*delta*k)*N.sinc(2*delta*l)*pi**3
            cosqr=N.cos(2*pi*1*(h*xi+l*zi));
            #fsum=(Aj*cosqr).sum()
            fsum=(fq*cosqr).sum()
            #print xi,zi,'sum',fsum 
            #for i=1:n
            #    currf=fq(i);
            #    h=Qs(i,1);%/q(i);
            #    k=Qs(i,2);%/q(i);
            #    l=Qs(i,3);%/q(i);
            #    Aj=currf*N.sinc(2*pi*delta*h)*sinc(2*pi*delta*k)*sinc(2*pi*delta*l);   
            #    eiqr=cos(2*pi*1*(h*x(xi)/aa+l*z(zi)/cc));
            #    fsum=fsum+Aj*eiqr;
            #    %fprintf('h=%f k=%f l=%f currf=%f Aj%f
            #    %eiqr=%f\n',h*q(i),k*q(i),l*q(i),currf,Aj,eiqr)
            #end     
            #P[zi,xi]=fsum
            P[xia,zia]=fsum
                #fprintf(fid,'%3.5g  %3.5g  %3.5g  \n',xi,zi,P(zi,xi));
    #print P
    print P.shape
    P_fou=dct.dct2(P)
    if 0:
        pylab.pcolor(X,Z,P_fou)
        pylab.show()
    return X,Z,P






def fourier_p(h,k,l,P,xstep=0.01,zstep=0.01):
    x=N.arange(0.0,1.0,xstep)
    z=N.arange(0.0,1.0,zstep)
    xn=len(x)
    zn=len(z)
    X,Z=N.meshgrid(x,z)
    fsum=0.0
    delta=.065;
    #%delta=.1;
    delta=.035;
    #outfile='c:\structfactors_density.dat'
    #fid=fopen([outfile],'w');
    fsum=0.0
    for xia in range(xn):
        for zia in range(zn):
            xi=x[xia]
            zi=z[zia]
            #Aj=fq*N.sinc(2*delta*h)*N.sinc(2*delta*k)*N.sinc(2*delta*l)*pi**3
            cosqr=N.cos(2*pi*1*(h*xi+l*zi));
            fsum=fsum+P[xia,zia]*cosqr  
            #print xi,zi,'sum',fsum,cosqr, P[xia,zia]
            #if abs(P[xia,zia])>0:
            #    print xi,zi,P[xia,zia]
    return 2*fsum/xn/zn

def transform_p(p,Mx,Mz,M):
    pup_n=len(p)/2
    pup=p[0:M]
    pdown=p[M::]
    #print 'pup',pup.shape
    #print 'pdown',pdown.shape
    p_up=pup.reshape(Mx,Mz)
    p_down=pdown.reshape(Mx,Mz)
    return p_up,p_down
     
     

def chisq(p,h,k,l,fq,fqerr,xstep=0.01,zstep=0.01):
    M=len(p)/2
    Mx=1.0/xstep
    Mz=1.0/zstep
    #print M,Mx,Mz
    fsum_up=N.zeros(h.shape)
    fsum_down=N.zeros(h.shape)
    chi=N.zeros(h.shape)
    P_up,P_down=transform_p(p,Mx,Mz,M)
    for i in range(len(h)):
        fsum_up[i]=fourier_p(h[i],k[i],l[i],P_up)
        fsum_down[i]=fourier_p(h[i],k[i],l[i],P_down)
        fmodel=fsum_up[i]-fsum_down[i]
        chi[i]=(fmodel-fq[i])**2/fqerr[i]**2
        #print h[i],k[i],l[i],fq[i],chi[i]
    
    return chi.sum()
    

def chisq_grad(p,h,k,l,fq,fqerr,xstep=0.01,zstep=0.01):
    M=len(p)/2
    Mx=1.0/xstep
    Mz=1.0/zstep
    #print M,Mx,Mz
    fsum_up=N.zeros(h.shape)
    fsum_down=N.zeros(h.shape)
    fmodel=N.zeros(h.shape)
    
    chi=N.zeros(h.shape)
    P_up,P_down=transform_p(p,Mx,Mz,M)

    x=N.arange(0.0,1.0,xstep)
    z=N.arange(0.0,1.0,zstep)
    xn=len(x)
    zn=len(z)
    grad=N.zeros(P_up.shape)
    for i in range(len(h)):
        fsum_up[i]=fourier_p(h[i],k[i],l[i],P_up)
        fsum_down[i]=fourier_p(h[i],k[i],l[i],P_down)
        fmodel[i]=fsum_up[i]-fsum_down[i]
    for xia in range(xn):
        for zia in range(zn):
            xi=x[xia]
            zi=z[zia]
            #Aj=fq*N.sinc(2*delta*h)*N.sinc(2*delta*k)*N.sinc(2*delta*l)*pi**3
            #cosqr=N.cos(2*pi*1*(h*xi+l*zi));
            #fsum=fsum+P[xia,zia]*cosqr  
            chi=0
            for i in range(len(h)):
                #fsum_up[i]=fourier_p(h[i],k[i],l[i],P_up)
                #fsum_down[i]=fourier_p(h[i],k[i],l[i],P_down)
                #fmodel=fsum_up[i]-fsum_down[i]
                cosqr=N.cos(2*pi*1*(h[i]*xi+l[i]*zi));
                chi=chi+2*(fmodel[i]-fq[i])/fqerr[i]**2*cosqr
                #print i
            grad[xia,zia]=chi.sum()
            
    grad_up=grad.flatten()
    grad_down=-grad_up
    grad=N.concatenate((grad_up,grad_down))
        #print h[i],k[i],l[i],fq[i],chi[i]
    
    return grad


def chisq_hessian(p,h,k,l,fq,fqerr,xstep=0.01,zstep=0.01):
    M=len(p)/2
    Mx=1.0/xstep
    Mz=1.0/zstep
    #print M,Mx,Mz
    fsum_up=N.zeros(h.shape)
    fsum_down=N.zeros(h.shape)
    fmodel=N.zeros(h.shape)
    
    chi=N.zeros(h.shape)
    P_up,P_down=transform_p(p,Mx,Mz,M)

    x=N.arange(0.0,1.0,xstep)
    z=N.arange(0.0,1.0,zstep)
    xn=len(x)
    zn=len(z)
    grad=N.zeros(P_up.shape)
    for i in range(len(h)):
        fsum_up[i]=fourier_p(h[i],k[i],l[i],P_up)
        fsum_down[i]=fourier_p(h[i],k[i],l[i],P_down)
        fmodel[i]=fsum_up[i]-fsum_down[i]
    for xia in range(xn):
        for zia in range(zn):
            xi=x[xia]
            zi=z[zia]
            #Aj=fq*N.sinc(2*delta*h)*N.sinc(2*delta*k)*N.sinc(2*delta*l)*pi**3
            #cosqr=N.cos(2*pi*1*(h*xi+l*zi));
            #fsum=fsum+P[xia,zia]*cosqr  
            chi=0
            for i in range(len(h)):
                #fsum_up[i]=fourier_p(h[i],k[i],l[i],P_up)
                #fsum_down[i]=fourier_p(h[i],k[i],l[i],P_down)
                #fmodel=fsum_up[i]-fsum_down[i]
                cosqr=N.cos(2*pi*1*(h[i]*xi+l[i]*zi));
                chi=chi+2*(fmodel[i]-fq[i])/fqerr[i]**2*cosqr
                #print i
            grad[xia,zia]=chi.sum()
            
    grad_up=grad.flatten()
    grad_down=-grad_up
    grad=N.concatenate((grad_up,grad_down))
        #print h[i],k[i],l[i],fq[i],chi[i]
    
    return grad





if __name__=="__main__":
    
    myfilestr=r'c:\structfactors.dat'
    bob=N.loadtxt(myfilestr)
    print bob
    print bob.shape
    h,k,l,fq,fqerr=N.loadtxt(myfilestr).T
    print h
    X,Z,P=plotdensity(h,k,l,fq)
    pu=P.flatten()
    pd=N.zeros(len(pu))
    p=N.concatenate((pu,pd))
    print len(pu)
    print len(pd)
    print len(p)
    chi=chisq(p,h,k,l,fq,fqerr)
    print chi
    grad=chisq_grad(p,h,k,l,fq,fqerr)
    print 'gradient',grad
    if 0:
        pylab.pcolor(X,Z,P)
        pylab.show()    
    fsum=fourier_p(h[0],k[0],l[0],P)
    fsum=N.zeros(h.shape)
    if 0:
        for i in range(len(h)):
            fsum[i]=fourier_p(h[i],k[i],l[i],P)
            print h[i],k[i],l[i],fq[i],fsum[i],fq[i]/fsum[i]
    
    if 0:
        q=N.sqrt(h**2+l**2)  #Note, 107 is off by a factor of 2!
        pylab.plot(q,fsum,'s')
        #pylab.twinx()
        pylab.plot(q,fq,'s')
        pylab.show()
    if 0:
        myfilestr=r'c:\structfactors_density.dat'
        x,y,z=N.loadtxt(myfilestr).T
        bob=N.loadtxt(myfilestr).T
        print bob.shape
        pylab.pcolormesh(x,y,z)
        pylab.show()