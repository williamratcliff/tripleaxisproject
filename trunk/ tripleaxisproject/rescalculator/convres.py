import numpy as N
import scipy
import scipy.integrate as integrate
import lattice_calculator
import rescalc
import copy
pi=N.pi
def ConvRes(sqw,pref,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=[7,0]):
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
    detM=MMxx*Myy*Mzz*Mww


    tqz=1./N.sqrt(Mzz); tqx=1/N.sqrt(MMxx)
    tqyy=1./N.sqrt(Myy); tqyx=-Mxy/Myy/N.sqrt(MMxx)
    tqww=1./N.sqrt(Mww); tqwy=-Myw/Mww/N.sqrt(Myy)
    tqwx=-(Mxw/Mww-Myw/Mww*Mxy/Myy)/N.sqrt(MMxx)

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
    testsqw=sqw(H,K,L,W,p)
    modes=testsqw.shape[0]
    print 'testsqw ',testsqw
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
    print 'prefactor_shape ',prefactor.shape
    if METHOD=='fixed':
       M=ACCURACY
       print 'M ',M
       step1=pi/(2*M[0]+1); step2=pi/(2*M[1]+1)
       dd1=N.linspace(-pi/2+step1/2,pi/2-step1/2,(2*M[0]+1))
       dd2=N.linspace(-pi/2+step2/2,pi/2-step2/2,(2*M[1]+1))
       print 'modes ',modes
       print 'points ',points
       print 'step1 ',step1
       print 'step2 ',step2
       print 'prefactor ',prefactor
       print 'dd1 ',dd1
       print 'dd2',dd2
       convs=N.zeros((modes,npts),'float64')
       conv=N.zeros((npts,),'float64')
       #cx,cy,cw=N.mgrid(dd1,dd1,dd1)
       a=-pi/2+step1/2
       b=pi/2-step1/2
       dpts=N.complex(0,2*M[0]+1)
       [cx,cy,cw]=N.mgrid[a:b:dpts,a:b:dpts,a:b:dpts]
       #tx=N.tan(N.reshape(cx,(N.size(cx),)))
       tx=N.tan(cx.flatten())
       ty=N.copy(tx)
       tw=N.copy(tx)
       #N.tan(N.reshape(cy,(N.size(cy),))); N.tan(N.reshape(cw,(N.size(cw),)))
       print 'tx ',tx.shape
       tz=N.tan(dd2)
       norm=(1+tx**2)*(1+ty**2)*(1+tw**2)*N.exp(-0.5*(tx**2+ty**2))*N.exp(-0.5*(tw**2))
       normz=N.exp(-0.5*(tz**2))*(1+tz**2)
       print 'tx ', tx.shape
       print 'ty ',ty.shape
       print 'tw ',tw.shape
       print 'tz ',tz.shape
       print 'norm ',norm
       print 'normz ',normz
       print 'cx ',cx.shape
       for iz in range(N.size(tz)):
           for i in range(npts):
              dQ1=tqx[i]*tx
              dQ2=tqyy[i]*ty+tqyx[i]*tx
              dW=tqwx[i]*tx+tqwy[i]*ty+tqww[i]*tw
              dQ4=tqz[i]*tz[iz]
              H1=H[i]+dQ1*xvec[0,i]+dQ2*yvec[0,i]+dQ4*zvec[0,i]
              K1=K[i]+dQ1*xvec[1,i]+dQ2*yvec[1,i]+dQ4*zvec[1,i]
              L1=L[i]+dQ1*xvec[2,i]+dQ2*yvec[2,i]+dQ4*zvec[2,i]
              W1=W[i]+dW;
              print 'dQ1 ',dQ1
              print 'dQ3 ',dQ2
              print 'dW ',dW
              print 'dQ4 ',dQ4
              print 'H1 ',H1
              print 'K1 ',K1
              print 'L1 ',L1
              print 'W1 ',W1
              print 'H1_shape ',H1.shape
              myint=sqw(H1,K1,L1,W1,p)
              for pt in range(15):
                print 'myint H %f K %f L %f W %f int %f'%(H1[pt],K1[pt],L1[pt],W1[pt],myint[0,pt])
              for j in range(modes):
                  add=myint[j,:]*norm*normz[iz]
                  convs[j,i]=convs[j,i]+add.sum()
                  print 'convs ',convs
              #print 'conv ',conv.shape
              conv[i]=(convs[:,i]*prefactor[:,i]).sum()
           print 'conv',conv
       factor=step1**3*step2/N.sqrt(detM)
       print 'factor ',factor
       conv=conv*factor
       if M[1]==0:
            conv=conv*0.79788
       if M[0]==0:
            conv=conv*0.79788**3

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


