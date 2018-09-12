import sys,os
import numpy as N
import pylab
import random
import matplotlib
#import scipy.sandbox.delaunay as D
import matplotlib.delaunay as D
#import numpy.core.ma as ma
import matplotlib.numerix.ma as ma
from matplotlib.ticker import NullFormatter, MultipleLocator,MaxNLocator
#from scipy.signal.signaltools import convolve2d

def plot_nodes(tri):
    for nodes in tri.triangle_nodes:
        D.fill(x[nodes],y[nodes],'b')
    pylab.show()

def plot_data(xa,ya,za,fig,nfig,colorflag=False,convolveflag=False):

    cmap = pylab.cm.jet
    cmap.set_bad('w', 1.0)
    myfilter=N.array([[0.1,0.2,0.1],[0.2,0.8,0.2],[0.1,0.2,0.1]],'d') /2.0
    #if convolveflag:
    #    zout=convolve2d(za,myfilter,mode='same') #to convolve, or not to convolve...
    #else:
    #    zout=za
    zout=za
    #zima = ma.masked_where(N.isnan(zout),zout)
    zima=za

    ax=fig.add_subplot(1,1,nfig)
    pc=ax.pcolormesh(xa,ya,zima,shading='interp',cmap=cmap)  # working good!
#    pc=ax.imshow(zima,interpolation='bilinear',cmap=cmap)
    
    pmin=zima.min()
    pmax=zima.max()
    #pmin=0
    #pmax=700
    #pc.set_clim(0.0,660.0)
    pc.set_clim(pmin,pmax)



    if colorflag:
        #g=pylab.colorbar(pc,ticks=N.arange(0,675,100))
        #g=pylab.colorbar(pc,ticks=N.arange(pmin,pmax,10000))
        g=pylab.colorbar(pc)
        #print g
        #g.ticks=None
        #gax.yaxis.set_major_locator(MultipleLocator(40))
        #g.ticks(N.array([0,20,40,60,80]))

    return ax,g




def prep_data2(xt,yt,zorigt):
#    Data=pylab.load(r'c:\resolution_stuff\1p4K.iexy')
    #Data=pylab.load(filename)
    #xt=Data[:,2]
    #yt=Data[:,3]
    #zorigt=Data[:,0]
    x=xt[:,zorigt>0.0]
    y=yt[:,zorigt>0.0]
    z=zorigt[:,zorigt>0.0]
#    zorig=ma.array(zorigt)
    print('reached')
    threshold=0.0;
#    print zorigt < threshold
#    print N.isnan(zorigt)
#    z = ma.masked_where(zorigt < threshold , zorigt)
    print('where masked ', z.shape)
#should be commented out--just for testing
##    x = pylab.randn(Nu)/aspect
##    y = pylab.randn(Nu)
##    z = pylab.rand(Nu)
##    print x.shape
##    print y.shape
    # Grid
    
    xi, yi = N.mgrid[-5:5:100j,-5:5:100j]
    xi,yi=N.mgrid[x.min():x.max():.001,y.min():y.max():.001]
    #zi=matplotlib.mlab.griddata(x,y,z,xi,yi)
    # triangulate data
    tri = D.Triangulation(x,y)
    #print 'before interpolator'
    # interpolate data
    interp = tri.nn_interpolator(z)
    #print 'interpolator reached'
    zi = interp(xi,yi)
    # or, all in one line
    #    zi = Triangulation(x,y).nn_interpolator(z)(xi,yi)
#    return x,y,z
    return xi,yi,zi




        
def get_tokenized_line(myfile,returnline=[''],delimiter=None):
    lineStr=myfile.readline()
    returnline[0]=lineStr.rstrip()
    strippedLine=lineStr.lower().rstrip()
    tokenized=strippedLine.split(delimiter)
    return tokenized


def read_specfile(filename,data,fields):
    myfile=open(filename,'r')
    print('inner', filename)
    #column_headers=myfile.readline()
    i=0
    #get scan number
    tokenized=get_tokenized_line(myfile)
    scan_num=int(tokenized[1])
    #read header lines
    while 1:
        tokenized=get_tokenized_line(myfile)
        #print tokenized
        if tokenized[0]=="#l":
            break
    if data=={}:
        print('empty',filename)
        fields=tokenized[1:]
        #print 'fields',fields
        for field in fields:
            data[field]=[]
    while 1:
        #print i
        tokenized=get_tokenized_line(myfile)
        if tokenized==[]:
            break
        if tokenized[0]=="#r":
            continue
        for curr in range(len(fields)):
            #if fields[curr]=='h':
            #    print 'h',filename
            if 1:
                data[fields[curr]].append(float(tokenized[curr]))
            if 0:
                data[fields[curr]].append(tokenized[curr])
                
        i=i+1
    myfile.close()
    return data,fields

if __name__=="__main__":
    mydirectory=r"C:\Temp\anl\Ratcliff\BFO_110.scans"
    #myfilebase=r"BFO_110.out.0112"
    start,stop=(112,125) #222
    start,stop=(126,137)  #h,l  #220
    #start,stop=(159,179) # h,k
    #start,stop=(267,279) # k,l
    data={}
    fields=[]
    myfilebase=r"BFO_110.out.0"
    for i in range(start,stop):
        myfilebasen=myfilebase+str(i)
        myfilestr=os.path.join(mydirectory,myfilebasen)
        print('outer',myfilestr)
        data,fields=read_specfile(myfilestr,data,fields)
        #print 'data',data['h']
    print('done')
    #print data
    if 0:
        pylab.plot(data['h'],data['cyber8c'],'s')
        pylab.show()
    if 1:
        if 1:
            x=N.array(data['h'])
        if 0:
            x=N.array(data['k'])       
        if 1:
            y=N.array(data['l'])
        if 0:
            y=N.array(data['k'])
            
        
        z=N.array(data['cyber8c'])
        xind=N.argsort(x)
        #x=x[xind]
        #y=y[xind]
        #z=z[xind]
        xa,ya,za=prep_data2(x,y,z)
        fig=pylab.figure(figsize=(8,8))
        ax,g=plot_data(xa,ya,za,fig,1,colorflag=True)
        pylab.show()
        
        