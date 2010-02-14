import numpy as np
from numpy import pi,cos, sqrt, sin,exp,conj
from numpy.linalg import norm
I=np.complex(0,1)
import copy
#from enthought.mayavi import mlab


astar=2*pi/5.511;
bstar=2*pi/5.511
cstar=2*pi/12.136

#def readfiles():
 #filedir='c:\srfeas\20081212\';
#filehead='mag35';
#fileend='.bt9';
#%scans=[12 16 18 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73];

#fileszero='00';


#fileszero2='0';
#%scans2=10:17;
#%scans=[scans1 scans2];
#%scans=scans1;
#files=[];
#numfiles=[2,3,5,6,7];

#numfiles=[4:9];
#for i=1:length(numfiles)
#file=[filedir filehead fileszero num2str(numfiles(i)) fileend];
#files=[files;file];
#end

#numfiles=[10:17];
#for i=1:length(numfiles)
#file=[filedir filehead fileszero2 num2str(numfiles(i)) fileend];
#files=[files;file];
#end

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
 s=np.vstack((s1,s2,s3,s4,s5,s6,s7,s8))
 return s
 

#@mlab.show
#def draw_struct():
    #fig=mlab.figure()    
    #r=gen_fe()
    #s=gen_spins()
    ##view along z-axis
    #x=r[:,0]
    #y=r[:,1]
    #z=r[:,2]
    #u=s[:,0]
    #v=s[:,1]
    #w=s[:,2]
    
    ##print x.shape
    ##print y.shape
    ##print z.shape
    #pts_as=mlab.points3d(x,y,z,color=(0,0,1),colormap='gist_rainbow',figure=fig,scale_factor=.1)
    #mlab.quiver3d(x, y, z,u,v,w, line_width=3, scale_factor=.3,figure=fig)
    #outline=mlab.outline(figure=fig,extent=[0,1,0,1,0,1])
    #mlab.orientation_axes(figure=fig,xlabel='a',ylabel='b',zlabel='c')
    #print 'done'

    
def mgnfacFesquared(x):

#x=s=sin(theta)/lambda=modq/4/pi
 y=( 0.0706*exp(-35.008*x**2)+0.3589*exp(-15.358*x**2)
   +0.5819*exp(-5.561*x**2)-0.0114)**2
 return y
    
def calcstructure(pfit,Qs,correction):
 #positions of magnetic atoms in unit cell
 Mn3pt=gen_fe()
 momentsKMn3p=gen_spins()
 n,m=Qs.shape
 QA=copy.deepcopy(Qs)
 QA[:,0]=QA[:,0]*astar
 QA[:,1]=QA[:,1]*bstar
 QA[:,2]=QA[:,2]*cstar
 modq=sqrt(QA[:,0]*QA[:,0]+QA[:,1]*QA[:,1]+QA[:,2]*QA[:,2])
 Qn=copy.deepcopy(QA)
 Qn[:,0]=QA[:,0]/modq;
 Qn[:,1]=QA[:,1]/modq;
 Qn[:,2]=QA[:,2]/modq;   
 #calculate form factor
 magnfacMn3p=sqrt(mgnfacFesquared(modq/4/pi));
 #calculate structure factor
 F1Mn3p=np.zeros((n,3));
 expt=exp(I*2*pi*np.dot(Qs,Mn3pt.T))
 F1Mn3p=np.dot(expt,momentsKMn3p)
 F1Mn3p[:,0]=F1Mn3p[:,0]*magnfacMn3p;
 F1Mn3p[:,1]=F1Mn3p[:,1]*magnfacMn3p;
 F1Mn3p[:,2]=F1Mn3p[:,2]*magnfacMn3p;
 #fmt=np.zeros((n,1))
 F1=F1Mn3p;
 #find fperpendicular from F-F.q*q hat
 Fperp1=np.zeros(F1.shape,'float64')
 Fperp1[:,0]=F1[:,0]-(Qn[:,0]*F1[:,0]+Qn[:,1]*F1[:,1]+Qn[:,2]*F1[:,2])*Qn[:,0]
 Fperp1[:,1]=F1[:,1]-(Qn[:,0]*F1[:,0]+Qn[:,1]*F1[:,1]+Qn[:,2]*F1[:,2])*Qn[:,1]
 Fperp1[:,2]=F1[:,2]-(Qn[:,0]*F1[:,0]+Qn[:,1]*F1[:,1]+Qn[:,2]*F1[:,2])*Qn[:,2]
 fm1=np.zeros(F1.shape,'float64')
 fm1=(Fperp1[:,0]*conj(Fperp1[:,0])+Fperp1[:,1]*conj(Fperp1[:,1])+Fperp1[:,2]*conj(Fperp1[:,2]))      
 fmt=fm1;
 fm=fmt[0:n]*pfit[0]**2*correction
 return fm



if __name__=="__main__":
 print "main"
 Qs=np.array([[1,0,1],
              [1,0,3]
              ],
              'float64')
 Qs=np.array([[  1,     0,     7],
              [  3,     0,     1],
              [  3,     0,     3],
              [  3,     0,     5],
              [  1,     0,     1],
              [  1,     0,     3],
              [  1,     0,     5],
              [  1,     0,     7],
              [  1,     0,     9],
              [  1,     0,    11],
              [  3,     0,     7],
              [  3,     0,     9],
              [  5,     0,     1],
              [  5,     0,     3]],
              'float64')
 pfit=np.array([1.0],'float64')  
 correction=np.ones(Qs.shape[0])
 fm=calcstructure(pfit,Qs,correction)
 print 'fm',fm
 #r=gen_fe()
 #draw_struct()
