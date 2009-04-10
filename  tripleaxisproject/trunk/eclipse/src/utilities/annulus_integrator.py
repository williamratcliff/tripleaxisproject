import readncnr3 as readncnr
import numpy as N
import scriptutil as SU
import re
import simple_combine
import copy
import os
import pylab
pi=N.pi
from matplotlib.mlab import griddata


def grid(x,y,z):
    xmesh_step=.02
    ymesh_step=.02
    xrange=N.linspace(x.min(),x.max(),37)
    yrange=N.linspace(y.min(),y.max(),68-43)
    print xrange
    print yrange
    print x
    xi,yi=N.mgrid[x.min():x.max():xmesh_step,y.min():y.max():ymesh_step]
    #blah
    # triangulate data
    #tri = D.Triangulation(N.copy(x),N.copy(y))
    #print 'before interpolator'
    ## interpolate data
    #interp = tri.nn_interpolator(z)
    #print 'interpolator reached'
    #zi = interp(xi,yi)
    print xi.shape
    print yi.shape
    zi = griddata(x,y,z,xi,yi)
    return xi,yi,zi

    
    
def readfiles(flist):
    mydatareader=readncnr.datareader()
    Qx=N.array([])
    Qy=N.array([])
    Qz=N.array([])
    Counts=N.array([])
    for currfile in flist:
        #print currfile
        mydata=mydatareader.readbuffer(currfile)
        #print mydata.data.keys()
        a=mydata.metadata['lattice']['a']
        b=mydata.metadata['lattice']['b']
        c=mydata.metadata['lattice']['c']
        Qx=N.concatenate((Qx,N.array(mydata.data['qx'])*2*pi/a))
        Qy=N.concatenate((Qy,N.array(mydata.data['qy'])*2*pi/b))
        Qz=N.concatenate((Qz,N.array(mydata.data['qz'])*2*pi/c))
        Counts=N.concatenate((Counts,N.array(mydata.data['counts'])))
    #xa,ya,za=prep_data2(Qx,Qy,Counts);
    #print 'xa',xa.min(),xa.max()
    #print 'qx',Qx.min(),Qx.max()
        #x,y,z=grid(Qx,Qz,Counts)
    return Qx,Qz,Counts


def findpeaks(qx,qz,q,counts):
    
    qlist=N.linspace(q.min(),q.max(),100)
    counts_out=[]
    for l in range(len(qlist)):
        qsum=0
        q_spaced=qlist[l]
        num_counted=0
        for i in range(len(q)):
                if l>0:
                    if q[i]<=q_spaced and q[i]>qlist[l-1]:
                        qsum=qsum+counts[i]
                        num_counted=num_counted+1
                else:
                    if q[i]<=q_spaced and q[i]>0:
                        qsum=qsum+counts[i]
                        num_counted=num_counted+1
        if num_counted>0:
            counts_out.append(qsum/num_counted)
        else:
            print 'qsum',qsum
            counts_out.append(0)
            
    return qlist,counts_out


if __name__=='__main__':
    myfilebase='SrFeA0'
    myend='bt9'
    mydirectory=r'C:\srfeas\SrFeAsNi\Ni0p08\2009-03-diffraction'
    myfilebaseglob=myfilebase+'*.'+myend
    print range(43,69)
    flist=[]
    for i in range(43,69):
        currfile=os.path.join(mydirectory,myfilebase+str(i)+r"."+myend)
        #print 'currfile',currfile
        flist.append(currfile)
    #flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    qx,qz,counts=readfiles(flist)
    x,y,z=grid(qx,qz,counts)
    print qx.shape, qz.shape, counts.shape
    q=N.sqrt(qx**2+qz**2)
    qout,counts_out=findpeaks(qx,qz,q,counts)
    if 0:
        pylab.plot(qout,counts_out,'s')
        pylab.show()
    if 1:
        #QX,QZ=N.meshgrid(qx,qz)
        pylab.contourf(x,y,z,15)#,cmap=pylab.cm.jet)
    
        #pylab.pcolor(qx,qz,counts)
        pylab.colorbar()
        pylab.show()
    