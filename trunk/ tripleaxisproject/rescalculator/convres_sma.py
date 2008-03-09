import numpy as N
import scipy
import scipy.integrate as integrate
import lattice_calculator
import rescalc
import copy

pi=N.pi



def approx1(t):
    y=(t*.5641896)/(.5 + t**2)
    return y

def approx2(t,u):
    y=(t*(1.410474 + u*.5641896))/(.75+(u*(3.+u)))
    return y

def approx3(t):
    y=(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+0.5642236*t))))\
     /(16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))))
    return y

def approx4(t,u):
    y=(t*(36183.31-u*(3321.99-u*(1540.787-u*(219.031-u*(35.7668-u*(1.320522-u*.56419))))))\
        /(32066.6-u*(24322.8-u*(9022.23-u*(2186.18-u*(364.219-u*(61.5704-u*(1.84144-u))))))))
    return y

def voigt(x,a):
    nx=N.size(x);
    if N.size(a)==1:
         a=N.ones((nx,),'float64')
    y=N.zeros((nx,),'complex')
    t=a+N.complex(0,1)*x; ax=N.abs(x); s=ax + a; u=t**2
    for idx in range(10):
        print '%d a %f ax %f x %f s %f u %f'%(idx,a[idx],ax[idx],x[idx],s[idx],u[idx])
    good=N.where(a==0)[0]
    if good.shape[0]!=0:
        y[good]=N.exp(-x[good]**2)
    good=N.where((a>=15) | (s>=15))[0]
    if type(good)==type((1,2)):
        good=good[0]
        print 'innergood ',good
    if good.shape[0]!=0:
        print 'good', good
        y[good]=approx1(t[good])
    good=N.where((s<15) & (a<15) & (a>=5.5))[0]
    print 'good2',good
    if good.shape[0]!=0:
        y[good]=approx2(t[good],u[good])
    good=N.where((s<15) & (s>=5.5) & (a<5.5))[0]
    if good.shape[0]!=0:
        y[good]=approx2(t[good],u[good])
    good=N.where((s<5.5) & (a<5.5) & (a>=0.75))[0]
    if good.shape[0]!=0:
        y[good]=approx3(t[good])
    good=N.where((s<5.5) & (a>=0.195*ax-0.176) & (a<0.75))[0]
    if good.shape[0]!=0:
        y[good]=approx3(t[good])
    good=N.where(~((s<5.5) & (a>=0.195*ax-0.176)) & (a<0.75))[0]
    if good.shape[0]!=0:
        y[good]= N.exp(u[good])-approx4(t[good],u[good])
    y=N.real(y)
    return y




def ConvResSMA(sqw,pref,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=[7,0]):
    """"Numerically calculate the convolution of a user-defined cross-section
    function with the Cooper-Nathans resolution function for a 3-axis neutron
    scattering experiment. See manual for details."""

    #Get standard orthonormnal coordinate system, cell parameters and reciprocal cell parameters...
    xvec=myrescal.lattice_calculator.x
    yvec=myrescal.lattice_calculator.y
    zvec=myrescal.lattice_calculator.z
    #print 'xvec ',xvec[0,0]
    #print 'xvec ',xvec[1,0]
    #print 'yvec ',yvec
    #print 'zvec ',zvec
    npts=H.shape[0]
    #Calculate the resolution matrix...
    R0,RMS=myrescal.ResMatS(H,K,L,W,setup);
    #print 'shape ', RMS.shape
    #print RMS[3,3,:]
    #print type(RMS)
    Mzz=N.array(RMS[3,3,:])
    Mww=N.array(RMS[2,2,:])
    Mxx=N.array(RMS[0,0,:])
    Mxy=N.array(RMS[0,1,:])
    Mxw=N.array(RMS[0,2,:])
    Myy=N.array(RMS[1,1,:])
    Myw=N.array(RMS[1,2,:])
    Mxx=Mxx-Mxw**2/Mww
    Mxy=Mxy-Mxw*Myw/Mww
    Myy=Myy-Myw**2/Mww
    MMxx=Mxx-Mxy**2/Myy
    GammaFactor=N.sqrt(Mww/2); OmegaFactorx=Mxw/N.sqrt(2*Mww); OmegaFactory=Myw/N.sqrt(2*Mww)
    #detM=MMxx*Myy*Mzz*Mww
    detxy=N.sqrt(Mxx*Myy-Mxy**2)
    detz=N.sqrt(Mzz)

    tqz=1./detz; tqy=N.sqrt(Mxx)/detxy;
    tqxx=1./N.sqrt(Mxx); tqxy=Mxy/N.sqrt(Mxx)/detxy;

##    print 'Mzz ',Mzz
##    print 'Mww ',Mww
##    print 'Mxx ',Mxx
##    print 'Mxy ',Mxy
##    print 'Mxw ',Mxw
##    print 'Myy ',Myy
##    print 'Myw ',Myw
##    print 'detM ', detM
##    print 'tqz ',tqz
##    print 'tqx ', tqx
##    print 'tqyy ',tqyy
##    print 'tqyx ',tqyx
##    print 'tqww ',tqww
##    print 'tqwy ',tqwy
##    print 'tqwx ',tqwx


    #================================================================================================
    myargs=(H,K,L,W,p)
    dispersion,testsqw,WL0=sqw(H,K,L,p)
    modes=testsqw.shape[0]
    try:
        (modes,points)=testsqw.shape
    except:
        (modes,)=testsqw.shape
        points=1
    if pref==None:
       prefactor=N.ones((modes,points),'float64')
       bgr=0
    else:
        prefactor,bgr=pref(H,K,L,W,myrescal,p)

#    %========================= fix  method ===========================================================
    if METHOD=='fixed':
       M=ACCURACY
       step1=pi/(2*M[0]+1); step2=pi/(2*M[1]+1)
       dd1=N.linspace(-pi/2+step1/2,pi/2-step1/2,(2*M[0]+1))
       dd2=N.linspace(-pi/2+step2/2,pi/2-step2/2,(2*M[1]+1))
       convs=N.zeros((modes,npts),'float64')
       conv=N.zeros((npts,),'float64')
       #cx,cy,cw=N.mgrid(dd1,dd1,dd1)
       a=-pi/2+step1/2
       b=pi/2-step1/2
       dpts=N.complex(0,2*M[0]+1)
       cx,cy=N.mgrid[a:b:dpts,a:b:dpts]
       #tx=N.tan(N.reshape(cx,(N.size(cx),)))
       tx=N.tan(cx.reshape(N.size(cx),1,order='F').copy().T).flatten()
       ty=N.tan(cy.reshape(N.size(cx),1,order='F').copy().T).flatten()
       tz=N.tan(dd2)
       norm=(1+tx**2)*(1+ty**2)*N.exp(-0.5*(tx**2+ty**2))
       normz=N.exp(-0.5*(tz**2))*(1+tz**2)
       #print 'norm shape ', norm.shape,normz.shape,detxy.shape,detz.shape
       for iz in range(N.size(tz)):
           for i in range(npts):
              dQ1=tqxx[i]*tx-tqxy[i]*ty
              dQ2=tqy[i]*ty
              dQ4=tqz[i]*tz[iz]
              H1=H[i]+dQ1*xvec[0,i]+dQ2*yvec[0,i]+dQ4*zvec[0,i]
              K1=K[i]+dQ1*xvec[1,i]+dQ2*yvec[1,i]+dQ4*zvec[1,i]
              L1=L[i]+dQ1*xvec[2,i]+dQ2*yvec[2,i]+dQ4*zvec[2,i]
              disp,myint,WL=sqw(H1,K1,L1,p)
              #for idx in range(15):
              #    print '%d H1 %f K1 %f L1 %f myint %f'%(idx,H1[idx],K1[idx],L1[idx],myint[0,idx])
              for j in range(modes):
                 Gamma=WL[j,:]*GammaFactor[i]
                 Omega=GammaFactor[i]*(disp[j,:]-W[i])+OmegaFactorx[i]*dQ1+OmegaFactory[i]*dQ2
                 #print 'Gamma Shape Omega', Gamma.shape, Omega.shape,disp.shape,W.shape
                 vog=voigt(Omega,Gamma)
                 for idx in range(15):
                    #print '%d %d Gamma %f disp %f W %f Omega %f'%(idx,j, Gamma[idx],disp[j,idx],W[0],\
                    #    Omega[idx])
                    #print '%d norm %f normz %f detxy %f detz %f'%(idx,norm[idx],normz[0],detxy[0],detz[0])
                    print '%d %d voigt %f'%(idx,j,vog[idx])
                 add=myint[j,:]*voigt(Omega,Gamma)*norm*normz[iz]/detxy[i]/detz[i]
                 convs[j,i]=convs[j,i]+add.sum()
                 #print 'convs ',convs
              #print 'conv ',conv.shape
              conv[i]=(convs[:,i]*prefactor[:,i]).sum()
              #print 'conv ',conv[i]
       factor=step1**2*step2
       #print 'factor ',factor
       conv=conv*factor
       if M[1]==0:
            conv=conv*0.79788
       if M[0]==0:
            conv=conv*0.79788**2

##    %========================= mc  method =============================================================
##
##    if strcmp(METHOD,'mc')
##       found=1;
##       if isempty(ACCURACY)
##          ACCURACY=10;
##       end;
##       M=ACCURACY;
##       convs=zeros(modes,len);
##       conv=zeros(1,len);
##       for i=1:len
##       for MonteCarlo=1:M
##          r=randn(4,1000)*pi-pi/2;
##          cx=r(1,:); cy=r(2,:); cz=r(3,:);cw=r(4,:);
##          tx=tan(cx(1:end)); ty=tan(cy(1:end)); tz=tan(cz(1:end)); tw=tan(cw(1:end));
##          norm=exp(-0.5*(tx.^2+ty.^2+tz.^2+tw.^2)).*(1+tx.^2).*(1+ty.^2).*(1+tz.^2).*(1+tw.^2);
##          dQ1=tqx(i)*tx;
##          dQ2=tqyy(i)*ty+tqyx(i)*tx;
##          dW=tqwx(i)*tx+tqwy(i)*ty+tqww(i)*tw;
##          dQ4=tqz(i)*tz;
##          H1=H(i)+dQ1*xvec(1,i)+dQ2*yvec(1,i)+dQ4*zvec(1,i);
##          K1=K(i)+dQ1*xvec(2,i)+dQ2*yvec(2,i)+dQ4*zvec(2,i);
##          L1=L(i)+dQ1*xvec(3,i)+dQ2*yvec(3,i)+dQ4*zvec(3,i);
##          W1=W(i)+dW;
##          int=feval(sqw,H1,K1,L1,W1,p);
##          for j=1:modes
##              add=int(j,:).*norm;
##              convs(j,i)=convs(j,i)+sum(add);
##          end
##          conv(i)=sum( convs(:,i).*prefactor(:,i) );
##      end;end;
##       conv=conv/M/1000*pi^4./sqrt(detM);
##    end;
##    %==================================================================================================
##
##    if found==0:
##       error('??? Error in ConvRes: Unknown convolution method! Valid options are: "fix",  "mc".');

    conv=conv*R0
    conv=conv+bgr
    return conv


