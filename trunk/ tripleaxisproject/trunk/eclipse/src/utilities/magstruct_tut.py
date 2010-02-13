import numpy as np
from numpy import pi,cos, sqrt, sin,exp
from enthought.mayavi import mlab


astar=2*pi/5.511;
bstar=2*pi/5.511
cstar=2*pi/12.136

def magstruct(pfit,Qm,Int,Interr,correction,setup):
 Ifit=calcstructure(pfit,Qm,correction,setup)
 chisq=(Ifit-Int)/Interr;
 #[th1 phi2 th2]; 
 return chisq

def gen_fe():
 #generate atoms in the unit cell
 r1=[.25, .25, .25]
 r2=[.25, .75, .25]
 r3=[.75, .25, .25]
 r4=[.75,  .75, .25]
 
 r5=[.25, .25, .75]
 r6=[.25, .75, .75]
 r7=[.75, .25, .75]
 r8=[.75,  .75, .75]
 r=np.vstack((r1,r2,r3,r4,r5,r6,r7,r8))
 return r

def gen_spins(th_in=41):
 th=np.radians(th_in)
 #define the spins of the atoms in the cell
 s1=np.array([cos(th), sin(th), 0],'float64'); s1=s1/norm(s1);
 s2=np.array([cos(th), sin(th), 0],'float64'); s2=s2/norm(s2);
 s3=-s1
 s4=-s1
 s5=-s1
 s6=-s1
 s7=s1
 s8=s1
 

@mlab.show
def draw_struct():
    fig=mlab.figure()    
    r=gen_fe()   
    #view along z-axis
    pts_as=mlab.points3d(r[:,0],r[:0,1],r[:,2],color=(1,0,0),colormap='gist_rainbow',figure=fig,scale_factor=.1)
    outline=mlab.outline(figure=fig,extent=[0,1,0,1,-1,0])
    mlab.orientation_axes(figure=fig,xlabel='a',ylabel='b',zlabel='c')
    print 'done'

#def calcstructure(pfit,Qs,correction,setup):
 ##positions of magnetic atoms in unit cell

#%Mn3pt=[r1;r2;r3;r4;r5;r6;r7;r8;r9;r10;r11;r12;r13;r14;r15;r16];
#Mn3pt=[r1;r2;r3;r4;r5;r6;r7;r8];
#%momentsKMn3p=[s1;s2;s3;s4;s5;s6;s7;s8;s9;s10;s11;s12;s13;s14;s15;s16];
#momentsKMn3p=[s1;s2;s3;s4;s5;s6;s7;s8];

#%Mn3pt=[r6;r8;r9;r11;r13;r14;r15;r16];
#%momentsKMn3p=[s6;s8;s9;s11;s13;s14;s15;s16];
 #[n m]=size(Qs);
#% Qs
  #Qs=repmat(Qs,1,1);
#%  correction=repmat(correction,2,1);    

    #%Qs(n+1:2*n,1)=Qs(1:n,2); % hkl->(khl)
    #%Qs(n+1:2*n,2)=Qs(1:n,1); % hkl->(khl)
    #%Qs(2*n+1:3*n,1)=-Qs(1:n,1); % hkl->(-hkl)
    #%Qs(3*n+1:4*n,2)=-Qs(1:n,2); % hkl->(h-kl)



    #QA=Qs;

    #QA(:,1)=QA(:,1)*astar;
    #QA(:,2)=QA(:,2)*bstar;
    #QA(:,3)=QA(:,3)*cstar;
    #modq=sqrt(QA(:,1).*QA(:,1)+QA(:,2).*QA(:,2)+QA(:,3).*QA(:,3));
    #Qn(:,1)=QA(:,1)./modq;
    #Qn(:,2)=QA(:,2)./modq;
    #Qn(:,3)=QA(:,3)./modq;    
    
    #magnfacMn3p=sqrt(mgnfacFesquared(modq/4/pi));
    #%magnfacMn3p=(mgnfacFesquared(modq/4/pi));
  

    #magnfacMn3p=repmat(magnfacMn3p,1,1);
    

    #%calculate structure factor
#F1Mn3p=zeros(1*n,3);
#expt=exp(i*2*pi*Qs*Mn3pt');
#F1Mn3p=expt*momentsKMn3p;

#F1Mn3p(:,1)=F1Mn3p(:,1).*magnfacMn3p;
#F1Mn3p(:,2)=F1Mn3p(:,2).*magnfacMn3p;
#F1Mn3p(:,3)=F1Mn3p(:,3).*magnfacMn3p;
#fmt=zeros(1*n,1);


#F1=F1Mn3p;



#%F1=s1*exp(2*pi*i*q*r1)+s2*exp(i*2*pi*q*r2)+s3*exp(2*pi*i*q*r3)+s4*exp(i*2*pi*q*r4);


#%find fperpendicular from F-F.q*q hat
    #Fperp1(:,1)=F1(:,1)-(Qn(:,1).*F1(:,1)+Qn(:,2).*F1(:,2)+Qn(:,3).*F1(:,3)).*Qn(:,1);
    #Fperp1(:,2)=F1(:,2)-(Qn(:,1).*F1(:,1)+Qn(:,2).*F1(:,2)+Qn(:,3).*F1(:,3)).*Qn(:,2);
    #Fperp1(:,3)=F1(:,3)-(Qn(:,1).*F1(:,1)+Qn(:,2).*F1(:,2)+Qn(:,3).*F1(:,3)).*Qn(:,3);
     
    #fm1=(Fperp1(:,1).*conj(Fperp1(:,1))+Fperp1(:,2).*conj(Fperp1(:,2))+Fperp1(:,3).*conj(Fperp1(:,3)));      
    #fmt=fmt+(fm1);

#%    size(correction)
#%    size(fmt)
#%fm=(fmt(1:n))*pfit(1)^2.*correction;

#[M1,M2,S1,S2,A1,A2]=SpecGoToicp(Qs(:,1),Qs(:,2),Qs(:,3),0,setup);
#theta=abs(S2)/2;
#x=0.0037;
#%%x=pfit(3);
#x=0.0;
#ext=(1+.001*x*norm(abs(fmt))*1.6377^3./sin(theta*2)).^(-0.5);%*pfit(1)^4*.0907;
#ext=ext';
#%ext=1;
#fm=(fmt(1:n))*pfit(1)^2.*correction.*ext;


#%    fm=((fmt(1:n)+fmt(n+1:2*n)+fmt(2*n+1:3*n)+fmt(3*n+1:4*n))/4)*pfit(1)^2.*correction;
#return

if __name__=="__main__":
 print "main"
 #r=gen_fe()
 draw_struct()
